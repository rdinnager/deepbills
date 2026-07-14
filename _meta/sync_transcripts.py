"""
Copy this project's Claude Code conversation transcripts into the project folder.

Claude Code stores transcripts as one JSONL per session under
  ~/.claude/projects/<encoded-working-dir>/<session-uuid>.jsonl
keyed by WORKING DIRECTORY, not project. Sessions run from the shared Projects
root sit intermingled with unrelated hub sessions (email triage, daily brief,
project-status) in one store.

INCLUSION SIGNAL: a session belongs to this project if it *accessed files under
the project's folder* (Read/Edit/Write/Bash on those paths — real work reads the
project's CLAUDE.md and touches its files many times). This is far more reliable
than matching the project *name* in prose: the daily-brief / email-triage /
project-status sessions merely mention every project, they don't do file work.
Keeping inclusion = file-access (not name-mention) is what stops private hub
transcripts from leaking onto the shared Drive.

SIDE TRANSCRIPTS (subagents / teammates / workflow agents) COUNT AS THE SESSION.
A session's own JSONL is only the *parent* conversation. Delegated work is written
to sibling stores:
  <store>/<session-uuid>/subagents/*.jsonl     one per subagent / teammate
  <store>/<session-uuid>/workflows/*.json      workflow journals
and those side dirs may be filed under a *different* store than the parent (a
subagent that ran with cwd=<project> lands in the project's own store even when
the parent ran from the root). A parent that delegates all its file work to
teammates shows almost no path hits of its own, so counting only the parent
transcript makes exactly the biggest sessions invisible. We therefore:
  - pool path-access hits over the parent AND every side transcript, and
  - copy the side transcripts alongside the parent (they are the actual record).

This script is PROJECT-AGNOSTIC: it infers the project name from its own location
(<project>/_meta/sync_transcripts.py), so every project ships a byte-identical
copy.

Rules:
  - dedicated working-dir stores (sessions launched from inside the project,
    incl. nested sub-folders): copy wholesale. Discovered by listing the store,
    not hand-computed, so relocated / nested encodings are picked up.
  - shared root store: include a session if it has >= THRESHOLD path-access
    references to <project>/ (pooled over parent + side transcripts; excludes a
    shallow one-file scan) and isn't a scheduled-task run.
  - ALLOWLIST / DENYLIST give manual override.
Idempotent: one file per session UUID, overwritten in place as the session grows
(no duplicate copies). Side transcripts go to transcripts/<uuid>/. Writes
transcripts/index.md.

Run: python sync_transcripts.py
"""
import os, re, json, shutil, glob

HOME = os.path.expanduser("~")
STORE = os.path.join(HOME, ".claude", "projects")

# --- self-locate: <project>/_meta/sync_transcripts.py -----------------------
META_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(META_DIR)
PROJECT = os.path.basename(PROJECT_DIR)
DEST = os.path.join(META_DIR, "transcripts")
os.makedirs(DEST, exist_ok=True)

# encoded working dir of the shared Projects root (constant for every project
# that lives directly under it). Claude Code encodes a path by replacing every
# non-alphanumeric char with '-'.
ROOT_ENC = "G--Shared-drives-COBL-Data-Projects"
ROOT_STORE = os.path.join(STORE, ROOT_ENC)


def enc(s):
    return re.sub(r"[^A-Za-z0-9]", "-", s)


# dedicated stores = any store dir whose encoded working dir is the project
# folder itself OR a folder nested inside it (e.g. Talks/Oz-2026 ->
# ...-Projects-Talks-Oz-2026). Discovered by listing, per the setup runbook.
DEDICATED_PREFIX = ROOT_ENC + "-" + enc(PROJECT)
DEDICATED_STORES = []
ALL_STORES = []
if os.path.isdir(STORE):
    for d in sorted(os.listdir(STORE)):
        p = os.path.join(STORE, d)
        if not os.path.isdir(p):
            continue
        ALL_STORES.append(p)
        if d == DEDICATED_PREFIX or d.startswith(DEDICATED_PREFIX + "-"):
            DEDICATED_STORES.append(p)
# EXTRA_STORES: manually add oddly-encoded / relocated stores that genuinely
# hold this project's sessions but don't match the prefix above.
EXTRA_STORES = []
DEDICATED_STORES += [s for s in EXTRA_STORES if os.path.isdir(s)]

# path-access references to the project folder (forward slash, or JSON-escaped
# Windows backslash) — i.e. a file operation inside <project>/, not a prose
# mention.
PATH_RE = re.compile(re.escape(PROJECT) + r"[\\/]")
THRESHOLD = 3        # min path-access lines to count as "did work in the project"
ALLOWLIST = set()    # force-include these session UUIDs
DENYLIST = set()     # force-exclude these session UUIDs


def side_dirs(uuid):
    """Every <store>/<uuid>/ artifact dir, across ALL stores (a session's
    subagents can be filed under a different store than its parent)."""
    return [os.path.join(s, uuid) for s in ALL_STORES
            if os.path.isdir(os.path.join(s, uuid))]


def side_transcripts(uuid):
    """Subagent/teammate JSONLs belonging to this session, across all stores."""
    out = []
    for d in side_dirs(uuid):
        out += sorted(glob.glob(os.path.join(d, "subagents", "*.jsonl")))
    return out


def path_hits(path, max_bytes=128 * 1024 * 1024):
    hits, read = 0, 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            read += len(line)
            if PATH_RE.search(line):
                hits += 1
            if read > max_bytes:
                break
    return hits


def session_hits(uuid, main_path):
    """Path-access hits pooled over the parent transcript and every subagent
    transcript — a parent that delegates its file work would otherwise look
    like it never touched the project."""
    h = path_hits(main_path)
    for sub in side_transcripts(uuid):
        h += path_hits(sub)
    return h


def is_scheduled(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            try:
                o = json.loads(line)
            except Exception:
                continue
            if o.get("type") == "user":
                c = o.get("message", {}).get("content")
                txt = c if isinstance(c, str) else (
                    next((p.get("text", "") for p in c if isinstance(p, dict)
                          and p.get("type") == "text"), "") if isinstance(c, list) else "")
                return txt.strip().startswith("<scheduled-task")
    return False


def summarize(path):
    """First prompt, line count, time span. Falls back to the `last-prompt`
    metadata record when a transcript carries no message lines (e.g. a session
    whose history was truncated), so it still shows up meaningfully."""
    first_prompt, n, first_ts, last_ts, meta_prompt = None, 0, None, None, None
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            n += 1
            try:
                o = json.loads(line)
            except Exception:
                continue
            ts = o.get("timestamp")
            if ts:
                first_ts = first_ts or ts
                last_ts = ts
            if o.get("type") == "last-prompt" and not meta_prompt:
                meta_prompt = o.get("lastPrompt")
            if first_prompt is None and o.get("type") == "user":
                c = o.get("message", {}).get("content")
                if isinstance(c, str):
                    first_prompt = c
                elif isinstance(c, list):
                    for part in c:
                        if isinstance(part, dict) and part.get("type") == "text":
                            first_prompt = part.get("text"); break
    if first_prompt is None and meta_prompt:
        first_prompt = "(no messages in transcript; last-prompt metadata) " + meta_prompt
    if first_prompt:
        first_prompt = " ".join(first_prompt.split())[:140]
    return first_prompt, n, first_ts, last_ts


def copy_side(uuid):
    """Mirror this session's subagent transcripts + workflow journals into
    transcripts/<uuid>/. Returns (n_files, total_bytes)."""
    n, size = 0, 0
    for d in side_dirs(uuid):
        for sub in ("subagents", "workflows"):
            src_dir = os.path.join(d, sub)
            if not os.path.isdir(src_dir):
                continue
            dst_dir = os.path.join(DEST, uuid, sub)
            os.makedirs(dst_dir, exist_ok=True)
            for src in sorted(glob.glob(os.path.join(src_dir, "*"))):
                if not os.path.isfile(src):
                    continue
                dst = os.path.join(dst_dir, os.path.basename(src))
                if (not os.path.exists(dst)) or os.stat(dst).st_size != os.stat(src).st_size:
                    shutil.copy2(src, dst)
                n += 1
                size += os.path.getsize(src)
    return n, size


def main():
    picked = {}   # uuid -> (src, hits)  hits None = dedicated store
    for ded in DEDICATED_STORES:
        for p in glob.glob(os.path.join(ded, "*.jsonl")):
            picked[os.path.basename(p)[:-6]] = (p, None)
    borderline = []
    if os.path.isdir(ROOT_STORE):
        for p in glob.glob(os.path.join(ROOT_STORE, "*.jsonl")):
            uuid = os.path.basename(p)[:-6]
            if uuid in DENYLIST or uuid in picked:
                continue
            h = session_hits(uuid, p)
            include = (uuid in ALLOWLIST) or (h >= THRESHOLD and not is_scheduled(p))
            if include:
                picked[uuid] = (p, h)
            elif h > 0:
                borderline.append((uuid, h))

    copied, skipped, rows = 0, 0, []
    for uuid, (src, hits) in sorted(picked.items()):
        dst = os.path.join(DEST, uuid + ".jsonl")
        if (not os.path.exists(dst)) or os.stat(dst).st_size != os.stat(src).st_size:
            shutil.copy2(src, dst); copied += 1
        else:
            skipped += 1
        n_side, side_bytes = copy_side(uuid)
        prompt, nlines, first_ts, last_ts = summarize(dst)
        rows.append((uuid, first_ts, last_ts, os.path.getsize(dst), nlines, prompt,
                     hits, n_side, side_bytes))

    n_dedicated = sum(1 for _, (_, h) in picked.items() if h is None)
    n_root = len(picked) - n_dedicated
    n_side_total = sum(r[7] for r in rows)
    rows.sort(key=lambda r: (r[1] or ""))
    with open(os.path.join(DEST, "index.md"), "w", encoding="utf-8") as f:
        f.write(f"# {PROJECT} — transcript index\n\n")
        f.write("Raw Claude Code session transcripts (JSONL), auto-collected by "
                f"`../sync_transcripts.py` (sessions that did file work under `{PROJECT}/`, "
                "counting work delegated to subagents).\n\n")
        f.write(f"{len(rows)} session(s).\n\n")
        for uuid, first_ts, last_ts, sz, nlines, prompt, hits, n_side, side_bytes in rows:
            tag = "dedicated-store" if hits is None else f"{hits} path-access refs"
            f.write(f"- **`{uuid}.jsonl`**  ({tag})  \n")
            f.write(f"  {(first_ts or '?')[:19]} → {(last_ts or '?')[:19]} · "
                    f"{sz/1e6:.1f} MB · {nlines} lines  \n")
            if n_side:
                f.write(f"  side transcripts: `{uuid}/` — {n_side} file(s), "
                        f"{side_bytes/1e6:.1f} MB (subagents + workflow journals)  \n")
            f.write(f"  _first prompt:_ {prompt or '(n/a)'}\n")
    print(f"[{PROJECT}] transcripts: {len(rows)} in project "
          f"({n_dedicated} dedicated-store, {n_root} root-store); "
          f"{copied} copied, {skipped} unchanged; {n_side_total} side transcript(s) "
          f"-> {DEST}")
    if borderline:
        print(f"[{PROJECT}] excluded (below threshold / prose-mention / shallow scan):")
        for uuid, h in sorted(borderline, key=lambda x: -x[1]):
            print(f"    {uuid}  ({h} path-access refs)")


if __name__ == "__main__":
    main()
