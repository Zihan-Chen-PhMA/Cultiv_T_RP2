from __future__ import annotations
import collections
import dataclasses
import heapq
import math
import time
from typing import Literal, cast, Any, AbstractSet, Union, Optional

import numpy as np
import pymatching
import sinter
import stim
from dem_parsor import *


class SO3Sampler2D(sinter.Sampler):
    def compiled_sampler_for_task(self, task: sinter.Task) -> CompiledSO3Sampler2D:
        return CompiledSO3Sampler2D.from_task(task)
    pass

class CompiledSO3Sampler2D(sinter.CompiledSampler):
    def __init__(self, dem: stim.DetectorErrorModel,
                 circuit: stim.Circuit,
                 post_select_det_ids: frozenset,
                 decoder: Optional[pymatching.Matching] = None):
        self.dem = dem
        self.num_detectors = dem.num_detectors
        self.circuit = circuit
        self.sampler = self.circuit.compile_detector_sampler()
        self.post_select_det_ids = post_select_det_ids
        self.decoder = decoder
        self._mask = np.array([i in self.post_select_det_ids 
                                for i in range(self.num_detectors)],dtype=np.uint8)
        if self.decoder is None:
            # self.decoder = pymatching.Matching(self.dem)
            raise ValueError('decoder unset')
        

    @staticmethod
    def from_task(task: sinter.Task) -> CompiledSO3Sampler2D:
        stim_circ = task.circuit
        dem_helper = DEM(stim_circ)
        dem = dem_helper.prune_post_selected()
        post_select_det_ids = frozenset(dem_helper.post_selection_mask)
        decoder = construct_decoder(task, d_rp=3)
        return CompiledSO3Sampler2D(dem,stim_circ,post_select_det_ids,decoder)
    

    def sample(self, suggested_shots) -> sinter.AnonTaskStats:
        t0 = time.monotonic()
        dets, obs = self.sampler.sample(shots=suggested_shots, separate_observables=True)
        keep_mask = ~np.any(dets & self._mask, axis=1)
        dets = dets[keep_mask]
        obs = np.reshape(obs[keep_mask],(-1))
        predicted_obs, soft_output_mono, soft_output = self.decoder.decode_batch_soft_output_2d(dets,return_weights=True)
        # predicted_obs, soft_output = self.decoder.decode_batch(dets,return_weights=True)
        soft_output_rounded = np.round(np.reshape(soft_output,(-1))*10/(2.303)).astype(np.int64)
        soft_output_mono_rounded =  np.round(np.reshape(soft_output_mono,(-1))*10/(2.303)).astype(np.int64)
        errs = np.reshape(predicted_obs,(-1)) ^ obs
        counter = collections.Counter()
        for gap_mono, gap, err in zip(soft_output_mono_rounded ,soft_output_rounded, errs):
            counter[f'E{gap_mono}|{gap}' if err else f'C{gap_mono}|{gap}'] += 1
        t1 = time.monotonic()

        return sinter.AnonTaskStats(
            shots=suggested_shots,
            errors=np.count_nonzero(errs),
            discards=suggested_shots - np.count_nonzero(keep_mask),
            seconds=t1 - t0,
            custom_counts=counter,
        )



class SO5Sampler2D(sinter.Sampler):
    def compiled_sampler_for_task(self, task: sinter.Task) -> CompiledSO5Sampler2D:
        return CompiledSO5Sampler2D.from_task(task)
    pass

class CompiledSO5Sampler2D(sinter.CompiledSampler):
    def __init__(self, dem: stim.DetectorErrorModel,
                 circuit: stim.Circuit,
                 post_select_det_ids: frozenset,
                 decoder: Optional[pymatching.Matching] = None):
        self.dem = dem
        self.num_detectors = dem.num_detectors
        self.circuit = circuit
        self.sampler = self.circuit.compile_detector_sampler()
        self.post_select_det_ids = post_select_det_ids
        self.decoder = decoder
        self._mask = np.array([i in self.post_select_det_ids 
                                for i in range(self.num_detectors)],dtype=np.uint8)
        if self.decoder is None:
            # self.decoder = pymatching.Matching(self.dem)
            raise ValueError('decoder unset')
        

    @staticmethod
    def from_task(task: sinter.Task) -> CompiledSO5Sampler2D:
        stim_circ = task.circuit
        dem_helper = DEM(stim_circ)
        dem = dem_helper.prune_post_selected()
        post_select_det_ids = frozenset(dem_helper.post_selection_mask)
        decoder = construct_decoder(task, d_rp=5)
        return CompiledSO5Sampler2D(dem,stim_circ,post_select_det_ids,decoder)
    

    def sample(self, suggested_shots) -> sinter.AnonTaskStats:
        t0 = time.monotonic()
        dets, obs = self.sampler.sample(shots=suggested_shots, separate_observables=True)
        keep_mask = ~np.any(dets & self._mask, axis=1)
        dets = dets[keep_mask]
        obs = np.reshape(obs[keep_mask],(-1))
        predicted_obs, soft_output_mono, soft_output = self.decoder.decode_batch_soft_output_2d(dets,return_weights=True)
        # predicted_obs, soft_output = self.decoder.decode_batch(dets,return_weights=True)
        soft_output_rounded = np.round(np.reshape(soft_output,(-1))*10/(2.303)).astype(np.int64)
        soft_output_mono_rounded =  np.round(np.reshape(soft_output_mono,(-1))*10/(2.303)).astype(np.int64)
        errs = np.reshape(predicted_obs,(-1)) ^ obs
        counter = collections.Counter()
        for gap_mono, gap, err in zip(soft_output_mono_rounded ,soft_output_rounded, errs):
            counter[f'E{gap_mono}|{gap}' if err else f'C{gap_mono}|{gap}'] += 1
        t1 = time.monotonic()

        return sinter.AnonTaskStats(
            shots=suggested_shots,
            errors=np.count_nonzero(errs),
            discards=suggested_shots - np.count_nonzero(keep_mask),
            seconds=t1 - t0,
            custom_counts=counter,
        )




def construct_decoder(task: sinter.Task, d_rp: int) -> pymatching.Matching:
    dem = DEM(task.circuit)

    pruned_dem = dem.prune_post_selected()

    unpost_det_ids = []
    unpost_time = []
    for i in range(dem.num_dets):
        if i not in dem.post_selection_mask:
            unpost_det_ids.append(i)
            unpost_time.append(round(dem.det_id_coords[i][-1]))

    growth_time = min(unpost_time)

    x_inner_circ_dets = []
    z_inner_circ_dets = []
    x_outer_circ_dets = []
    z_outer_circ_dets = []

    x_upper_bd_dets = []
    x_lower_bd_dets = []
    z_left_bd_dets = []
    z_right_bd_dets = []



    for det_id in unpost_det_ids:
        det_id_coords = dem.det_id_coords[det_id]
        if det_id_coords[-1] == growth_time:
            # x on horizontal line (diagonal tile excluded)
            if np.abs(det_id_coords[1]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[0]) < d_rp//2 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 0: 
                x_inner_circ_dets.append(det_id)
                # print(det_id_coords,'hori')
            # x on vertical line (diagonal tile excluded)
            elif np.abs(det_id_coords[0]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[1]) < d_rp//2 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 0:
                x_inner_circ_dets.append(det_id)
                # print(det_id_coords,'vert')
            # x on diagonal line
            elif np.abs(det_id_coords[0]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[1]) == d_rp//2 + 0.5 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 0:
                x_inner_circ_dets.append(det_id)
                # print(det_id_coords,'diag')
            # z on horizontal line (diagonal tile excluded)
            elif np.abs(det_id_coords[1]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[0]) < d_rp//2 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 1:
                z_inner_circ_dets.append(det_id)
                # print(det_id_coords,'zhori')
            # z on vertical line (diagonal tile excluded) 
            elif np.abs(det_id_coords[0]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[1]) < d_rp//2 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 1:
                z_inner_circ_dets.append(det_id)
                # print(det_id_coords,'zvert')
            elif np.abs(det_id_coords[0]) == d_rp//2 + 0.5 \
                and np.abs(det_id_coords[1]) == d_rp//2 + 0.5 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 1:
                z_inner_circ_dets.append(det_id)
                # print(det_id_coords,'zdiag')
            # x on vertical line 
            elif np.abs(det_id_coords[0]) == d_rp//2 + 1.5 \
                and np.abs(det_id_coords[1]) < d_rp//2 + 1 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 0:
                x_outer_circ_dets.append(det_id)
                # print(det_id_coords,'xout')
            # z on horizontal line (outer)
            elif np.abs(det_id_coords[1]) == d_rp//2 + 1.5 \
                and np.abs(det_id_coords[0]) < d_rp//2 + 1 \
                and (det_id_coords[0]-det_id_coords[1]) % 2 == 1:
                z_outer_circ_dets.append(det_id)
                # print(det_id_coords,'zout')


    coords_to_det_id = {}
    for det_id in unpost_det_ids:
        det_id_coords = dem.det_id_coords[det_id]
        coords_to_det_id[tuple(det_id_coords)] = det_id
        try:
            err_dets_list = dem.det_id_connection[det_id]
        except:
            continue
        for err_dets in err_dets_list:
            if len(err_dets) == 1:
                if (det_id_coords[0] - det_id_coords[1]) % 2 == 0:
                    if det_id_coords[1] < 0:
                        x_upper_bd_dets.append(det_id)
                        # print(det_id_coords,'xupper')
                        break
                    elif det_id_coords[1] > 0:
                        x_lower_bd_dets.append(det_id)
                        # print(det_id_coords,'xlower')
                        break
                if (det_id_coords[0] - det_id_coords[1]) % 2 == 1:
                    if det_id_coords[0] < 0:
                        z_left_bd_dets.append(det_id)
                        # print(det_id_coords,'zleft')
                        break
                    elif det_id_coords[0] > 0:
                        z_right_bd_dets.append(det_id)
                        # print(det_id_coords,'zright')
                        break

    x_upper_index = dem.num_dets
    x_lower_index = x_upper_index + 1
    z_left_index = x_lower_index + 1
    z_right_index = z_left_index + 1

    new_index_temp = z_right_index + 1

    decoder = pymatching.Matching(pruned_dem)
    decoder.SO_calculator_setup()
    for bd_index in [x_upper_index,x_lower_index,
                    z_left_index,z_right_index]:
        decoder.add_boundary_node_SO(bd_index)


    cyc_err_coords_set: set[frozenset[tuple[float]]] = set()

    for det_id in z_inner_circ_dets + z_outer_circ_dets:
        det_coords = tuple(dem.det_id_coords[det_id])
        err_dets_list = dem.det_id_connection[det_id]
        for err_dets in err_dets_list:
            if len(err_dets) == 1:
                continue
            err_det = None
            for id in err_dets:
                if id != det_id:
                    err_det = id
            err_coords = tuple(dem.det_id_coords[err_det])
            if err_det in z_inner_circ_dets + z_outer_circ_dets:
                if distance(det_coords,err_coords) \
                    < distance(det_coords,inv(err_coords)):
                    cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
                    cyc_err_coords_set.add(frozenset([inv(det_coords),
                                                    inv(err_coords)]))
                elif distance(det_coords,err_coords) \
                    > distance(det_coords,inv(err_coords)):
                    cyc_err_coords_set.add(frozenset([det_coords,inv(err_coords)]))
                    cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
                else:
                    raise ValueError('strange')
            elif err_coords[0] == det_coords[0] \
                and err_coords[-1] == det_coords[-1]+1 \
                and (err_coords[1]-det_coords[1]) == 2: # hook err
                cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
            elif err_coords[0] == inv(det_coords)[0] \
                and err_coords[-1] == inv(det_coords)[-1]+1 \
                and (err_coords[1]-inv(det_coords)[1]) == 2:
                cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
            else:
                if distance(det_coords,err_coords) < distance(inv(det_coords),err_coords):
                    cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
                elif distance(det_coords,err_coords) > distance(inv(det_coords),err_coords):
                    cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
                else:
                    raise ValueError('strange')



    for det_id in x_inner_circ_dets + x_outer_circ_dets:
        det_coords = tuple(dem.det_id_coords[det_id])
        err_dets_list = dem.det_id_connection[det_id]
        for err_dets in err_dets_list:
            if len(err_dets) == 1:
                continue
            err_det = None
            for id in err_dets:
                if id != det_id:
                    err_det = id
            err_coords = tuple(dem.det_id_coords[err_det])
            if err_det in x_inner_circ_dets + x_outer_circ_dets:
                if distance(det_coords,err_coords) \
                    < distance(det_coords,inv(err_coords)):
                    cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
                    cyc_err_coords_set.add(frozenset([inv(det_coords),
                                                    inv(err_coords)]))
                elif distance(det_coords,err_coords) \
                    > distance(det_coords,inv(err_coords)):
                    cyc_err_coords_set.add(frozenset([det_coords,inv(err_coords)]))
                    cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
                else:
                    raise ValueError('strange')
            elif err_coords[1] == det_coords[1] \
                and err_coords[-1] == det_coords[-1]+1 \
                and (err_coords[0]-det_coords[0]) == 2: # hook err
                cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
            elif err_coords[1] == inv(det_coords)[1] \
                and err_coords[-1] == inv(det_coords)[-1]+1 \
                and (err_coords[0]-inv(det_coords)[0]) == 2:
                cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
            else:
                if distance(det_coords,err_coords) < distance(inv(det_coords),err_coords):
                    cyc_err_coords_set.add(frozenset([det_coords,err_coords]))
                elif distance(det_coords,err_coords) > distance(inv(det_coords),err_coords):
                    cyc_err_coords_set.add(frozenset([inv(det_coords),err_coords]))
                else:
                    raise ValueError('strange')


    cyc_coords_to_det = {}
    cyc_image_coords_to_det = {}
    cyc_real_image_pair = set()
    cyc_real2image = {}
    cyc_image2real = {}

    for det_id in z_inner_circ_dets + z_outer_circ_dets \
                + x_inner_circ_dets + x_outer_circ_dets:
        decoder.add_image_node_SO(det_id,new_index_temp)
        coords = tuple(dem.det_id_coords[det_id])
        cyc_coords_to_det[coords] = det_id
        cyc_image_coords_to_det[inv(coords)] = new_index_temp
        cyc_real_image_pair.add(frozenset([det_id,new_index_temp]))
        cyc_real2image[det_id] = new_index_temp
        cyc_image2real[new_index_temp] = det_id
        new_index_temp += 1

    cyc_real_coords = set(cyc_coords_to_det.keys())
    cyc_image_coords = set(cyc_image_coords_to_det.keys())

    for bd_det_id in z_left_bd_dets:
        decoder.add_boundary_edge_SO(bd_det_id,z_left_index)
        if bd_det_id in z_inner_circ_dets + z_outer_circ_dets:
            decoder.add_boundary_edge_to_image(bd_det_id,cyc_real2image[bd_det_id],
                                            z_right_index)
    for bd_det_id in z_right_bd_dets:
        decoder.add_boundary_edge_SO(bd_det_id,z_right_index)
        if bd_det_id in z_inner_circ_dets + z_outer_circ_dets:
            decoder.add_boundary_edge_to_image(bd_det_id,cyc_real2image[bd_det_id],
                                            z_left_index)
    for bd_det_id in x_lower_bd_dets:
        decoder.add_boundary_edge_SO(bd_det_id,x_lower_index)
        if bd_det_id in x_inner_circ_dets + x_outer_circ_dets:
            decoder.add_boundary_edge_to_image(bd_det_id,cyc_real2image[bd_det_id],
                                            x_upper_index)
    for bd_det_id in x_upper_bd_dets:
        decoder.add_boundary_edge_SO(bd_det_id,x_upper_index)
        if bd_det_id in x_inner_circ_dets + x_outer_circ_dets:
            decoder.add_boundary_edge_to_image(bd_det_id,cyc_real2image[bd_det_id],
                                            x_lower_index)

    while cyc_err_coords_set:
        cyc_err_coords = cyc_err_coords_set.pop()
        cyc_err_coords = list(cyc_err_coords)
        # print(cyc_err_coords)
        coord_1 = cyc_err_coords[0]
        coord_2 = cyc_err_coords[1]
        if coord_1 in cyc_real_coords:
            if coord_2 in cyc_real_coords:
                decoder.copy_edge_to_image(cyc_coords_to_det[coord_1],
                                        cyc_coords_to_det[coord_2],
                            cyc_image_coords_to_det[inv(coord_1)],
                            cyc_image_coords_to_det[inv(coord_2)])
                cyc_err_coords_set.remove(frozenset([
                        inv(coord_1), inv(coord_2)]))
            elif coord_2 in cyc_image_coords:
                decoder.redirect_edge_to_image(cyc_coords_to_det[inv(coord_2)],
                                cyc_image_coords_to_det[coord_2],
                                cyc_coords_to_det[coord_1])
                decoder.copy_edge_to_image(cyc_coords_to_det[coord_1],
                                        cyc_image_coords_to_det[coord_2],
                                cyc_image_coords_to_det[inv(coord_1)],
                                cyc_coords_to_det[inv(coord_2)])
                cyc_err_coords_set.remove(frozenset([inv(coord_1),inv(coord_2)]))
            else:
                pass
        elif coord_1 in cyc_image_coords:
            if coord_2 in cyc_real_coords:
                decoder.redirect_edge_to_image(cyc_coords_to_det[inv(coord_1)],
                                    cyc_image_coords_to_det[coord_1],
                                    cyc_coords_to_det[coord_2])
                decoder.copy_edge_to_image(cyc_image_coords_to_det[coord_1],
                                        cyc_coords_to_det[coord_2],
                                cyc_coords_to_det[inv(coord_1)],
                                cyc_image_coords_to_det[inv(coord_2)])
                cyc_err_coords_set.remove(frozenset([inv(coord_1),inv(coord_2)]))
            elif coord_2 in cyc_image_coords:
                decoder.copy_edge_to_image(cyc_coords_to_det[inv(coord_1)],
                                        cyc_coords_to_det[inv(coord_2)],
                                    cyc_image_coords_to_det[coord_1],
                                    cyc_image_coords_to_det[coord_2])
                cyc_err_coords_set.remove(frozenset([inv(coord_1),inv(coord_2)]))
            else:
                decoder.redirect_edge_to_image(cyc_coords_to_det[inv(coord_1)],
                                        cyc_image_coords_to_det[coord_1],
                                        coords_to_det_id[coord_2])
        else:
            if coord_2 in cyc_real_coords:
                pass
            else:
                decoder.redirect_edge_to_image(cyc_coords_to_det[inv(coord_2)],
                                        cyc_image_coords_to_det[coord_2],
                                        coords_to_det_id[coord_1])

    decoder.add_cycle_endpoints_pair_SO(z_left_index,z_right_index)
    decoder.add_cycle_endpoints_pair_SO(x_lower_index,x_upper_index)
    for det_id in cyc_real2image:
        if det_id in z_inner_circ_dets + x_inner_circ_dets:
            decoder.add_cycle_endpoints_pair_mono_SO(det_id,cyc_real2image[det_id])


    return decoder            