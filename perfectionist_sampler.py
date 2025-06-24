import time
import collections

import numpy as np
import sinter


class PerfectionistSampler(sinter.Sampler):
    """Predicts obs aren't flipped. Discards shots with any detection events."""
    def compiled_sampler_for_task(self, task: sinter.Task) -> sinter.CompiledSampler:
        return CompiledPerfectionistSampler(task)


class CompiledPerfectionistSampler(sinter.CompiledSampler):
    def __init__(self, task: sinter.Task):
        self.stim_sampler = task.circuit.compile_detector_sampler()
        self.dem = task.circuit.detector_error_model()
        self.num_detectors = self.dem.num_detectors
        self.id2coords = self.dem.get_detector_coordinates()
        self.stage2id: dict[int,list[int]] = {}
        self.stage2mask: dict[int,np.ndarray] = {} 
        for id in self.id2coords:
            coords = self.id2coords[id]
            if len(coords) != 4:
                raise NotImplementedError('all coords should have four components')
            try:
                self.stage2id[int(coords[-1])].append(id)
            except:
                self.stage2id[int(coords[-1])] = [id]
        self.stages = sorted(self.stage2id.keys())
        for stage in self.stages:
            ids = self.stage2id[stage]
            mask = np.array([id in ids for id in range(self.num_detectors)],dtype=np.uint8)
            self.stage2mask[stage] = mask
        

        

    def sample(self, max_shots: int) -> sinter.AnonTaskStats:
        t0 = time.monotonic()
        dets, obs = self.stim_sampler.sample(
            shots=max_shots,
            bit_packed=False,
            separate_observables=True,
        )
        num_shots = dets.shape[0]
        num_discards = 0
        counter = collections.Counter()
        for stage in self.stages:
            keep_mask = ~np.any(dets & self.stage2mask[stage], axis=1)
            counter[f'D{stage}'] = np.count_nonzero(keep_mask==0)
            num_discards += counter[f'D{stage}']
            dets = dets[keep_mask]
            obs = obs[keep_mask]

        errors = np.any(obs, axis=1)
        num_errors = np.count_nonzero(errors)
        t1 = time.monotonic()

        return sinter.AnonTaskStats(
            shots=num_shots,
            errors=num_errors,
            discards=num_discards,
            seconds=t1 - t0,
            custom_counts=counter,
        )
