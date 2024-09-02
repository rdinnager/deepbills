from model.sdf_net import SDFNet, LATENT_CODE_SIZE, LATENT_CODES_FILENAME
from util import device, standard_normal_distribution, ensure_directory
from dataset import dataset
import scipy
import numpy as np
from rendering import MeshRenderer
import time
import torch
from tqdm import tqdm
import cv2
import random
import sys

SAMPLE_COUNT = 30 # Number of distinct objects to generate and interpolate between
TRANSITION_FRAMES = 60

ROTATE_MODEL = False
USE_HYBRID_GAN = False
CATEGORY = None # Limit category when using the DeepSDF autodecoder, set to None to use all categories

SURFACE_LEVEL = 0.07 # 0.048 if USE_HYBRID_GAN else 0.011

sdf_net = SDFNet()
if USE_HYBRID_GAN:
    sdf_net.filename = 'hybrid_gan_generator.to'
sdf_net.load()
sdf_net.eval()

if USE_HYBRID_GAN:
    codes = standard_normal_distribution.sample((SAMPLE_COUNT + 1, LATENT_CODE_SIZE)).numpy()
else:
    latent_codes = torch.load(LATENT_CODES_FILENAME).detach().cpu().numpy()
    if CATEGORY is not None:
        dataset.load_labels()    
        indices = (dataset.labels == CATEGORY).nonzero()
        indices = random.sample(list(indices.cpu().numpy()), SAMPLE_COUNT + 1)
    else:
        indices = random.sample(list(range(latent_codes.shape[0])), SAMPLE_COUNT + 1)
    codes = latent_codes[indices, :]

codes[0, :] = codes[-1, :] # Make animation periodic
spline = scipy.interpolate.CubicSpline(np.arange(SAMPLE_COUNT + 1), codes, axis=0, bc_type='periodic')

def create_image_sequence():
    ensure_directory('images')
    frame_index = 0
    viewer = MeshRenderer(size=1080, start_thread=False)
    if CATEGORY is not None:
        viewer.model_color = dataset.get_color(CATEGORY)
    progress_bar = tqdm(total=SAMPLE_COUNT * TRANSITION_FRAMES)

    for sample_index in range(SAMPLE_COUNT):
        for step in range(TRANSITION_FRAMES):
            code = torch.tensor(spline(float(sample_index) + step / TRANSITION_FRAMES), dtype=torch.float32, device=device)
            if ROTATE_MODEL:
                viewer.rotation = (147 + frame_index / (SAMPLE_COUNT * TRANSITION_FRAMES) * 360 * 6, 40)
            viewer.set_mesh(sdf_net.get_mesh(code, voxel_resolution=128, sphere_only=True, level=SURFACE_LEVEL))
            image = viewer.get_image(flip_red_blue=True)
            cv2.imwrite("images/frame-{:05d}.png".format(frame_index), image)
            frame_index += 1
            progress_bar.update()
    
    print("\n\nUse this command to create a video:\n")
    print('ffmpeg -framerate 30 -i images/frame-%05d.png -c:v libx264 -profile:v high -crf 19 -pix_fmt yuv420p video.mp4')

def show_models():
    TRANSITION_TIME = 2
    viewer = MeshRenderer()

    if CATEGORY is not None:
        viewer.model_color = dataset.get_color(CATEGORY)

    while True:
        for sample_index in range(SAMPLE_COUNT):
            try:
                start = time.perf_counter()
                end = start + TRANSITION_TIME
                while time.perf_counter() < end:
                    progress = min((time.perf_counter() - start) / TRANSITION_TIME, 1.0)
                    if ROTATE_MODEL:
                        viewer.rotation = (147 + (sample_index + progress) / SAMPLE_COUNT * 360 * 6, 40)
                    code = torch.tensor(spline(float(sample_index) + progress), dtype=torch.float32, device=device)
                    viewer.set_mesh(sdf_net.get_mesh(code, voxel_resolution=32, sphere_only=True, level=SURFACE_LEVEL))
                
            except KeyboardInterrupt:
                viewer.stop()
                return

if 'save' in sys.argv:
    create_image_sequence()
else:
    show_models()
