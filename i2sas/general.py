import numpy as np


def get_volume_fraction(segmentation):
    particles_mask = ((segmentation > 0.0) | (np.bitwise_not(np.isnan(segmentation)))).astype(np.uint8)
    return np.sum(particles_mask) / np.prod(segmentation.shape)


def get_areas(segmentation):
    areas = []
    for inst in np.unique(segmentation):
        if inst == 0 or not np.isnan(inst):
            inst_mask = (segmentation == inst).astype(np.uint8)
            areas.append(np.sum(inst_mask))
    return np.array(areas)