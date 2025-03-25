import stim
from typing import Tuple, Union
import numpy as np
import copy


# process the stim.dem 


class Err_Unit():

    def __init__(self, 
                 err_location: stim.CircuitErrorLocation) -> None:
        '''
        Each err_unit is a specific error location


        flipped pauli product list[str]
        target qubit circuit_id list[int]
        instruction line: int
        tick offset: int
        error source: str
        error source probability: float
        '''
        self.err_location = err_location
        self.error_type, self.error_supp = self._pauli_product_extraction()
        self.error_source, self.error_source_probability = \
                                            self._error_source()
        self.instruction_line, self.instruction_index_list = \
                                            self._instruction()
        self.error_probability = self._err_probability()
        pass

    def _pauli_product_extraction(self) -> Tuple[
                                            list[str], list[int]]:
        flipped_pauli_product_list = self.err_location\
                                            .flipped_pauli_product
        pauli_type_list: list[str] = []
        pauli_supp_list: list[int] = []
        for gate_target in flipped_pauli_product_list:
            pauli_type_list.append(gate_target.gate_target.pauli_type)
            pauli_supp_list.append(gate_target.gate_target.qubit_value)
        return pauli_type_list, pauli_supp_list
    
    def _error_source(self) -> Tuple[str, float]:
        instruction_targets = self.err_location.instruction_targets
        error_source = instruction_targets.gate
        error_source_probability = instruction_targets.args[0]
        return error_source, error_source_probability
    
    def _instruction(self) -> Tuple[int, list[int]]:
        instruction_targets = self.err_location.instruction_targets
        instruction_index_list = [i for i in range(
            instruction_targets.target_range_start + 1,
            instruction_targets.target_range_end + 1)]
        stack_frame = self.err_location.stack_frames[0]
        instruction_line = stack_frame.instruction_offset + 1
        return instruction_line, instruction_index_list
    
    def _err_probability(self) -> float:
        if self.error_source == 'X_ERROR':
            return self.error_source_probability
        elif self.error_source == 'Z_ERROR':
            return self.error_source_probability
        elif self.error_source == 'DEPOLARIZE1':
            return self.error_source_probability/3
        elif self.error_source == 'DEPOLARIZE2':
            return self.error_source_probability/15
        else:
            raise ValueError('NO such err source')




    

class Fault_Unit():

    def __init__(self, fault_id: int, 
                 err_location_list: list[stim.CircuitErrorLocation],
                 syndrome_logic_list: list[stim.DemTarget]
                 ) -> None:
        '''
        each Fault_Unit is a collection of equivalent Err_Units (error locations)
        '''
        self.fault_id = fault_id
        self.err_unit_list = [Err_Unit(err_location) for err_location in 
                              err_location_list]
        self.err_probability_list = [err_unit.error_probability 
                                     for err_unit in self.err_unit_list]
        self.err_probability = self._err_probability()
        self.syndrome_logic_list = syndrome_logic_list
        self.syndrome_flag, self.flipped_detector_id_list = \
                                                    self._syndrome()
        self.logical_flag, self.flipped_logical_id_list = \
                                                    self._logical()
        self.err_type, self.err_terms, self.err_source = self._err_type()
        self.err_type_all = self._err_type_all()
        pass

    def _syndrome(self) -> Tuple[bool, list[int]]:
        detector_id_list = []
        for dem_target in self.syndrome_logic_list:
            if dem_target.is_relative_detector_id() == True:
                detector_id_list.append(dem_target.val)
        syndrome_flag = True
        if len(detector_id_list) == 0:
            syndrome_flag = False
        return syndrome_flag, detector_id_list

    def _logical(self) -> Tuple[bool, list[int]]:
        logic_id_list = []
        for dem_target in self.syndrome_logic_list:
            if dem_target.is_logical_observable_id() == True:
                logic_id_list.append(dem_target.val)
        logic_flag = False
        if len(logic_id_list) != 0:
            logic_flag = True
        return logic_flag, logic_id_list
    
    def _err_probability(self) -> float:
        prob_temp = 0
        for index, item in zip(range(len(self.err_probability_list)),
                               self.err_probability_list):
            err_prob_arr_copy = 1 - np.array(copy.deepcopy(
                                        self.err_probability_list))
            err_prob_arr_copy[index] = item
            prob_temp = prob_temp + np.prod(err_prob_arr_copy)
        return prob_temp

    def _err_type(self) -> Tuple[str,
                                 Union[list[str],None],
                                 Union[str,None]]:
        for err_unit in self.err_unit_list:
            if 'X' in err_unit.error_type and \
               'Y' not in err_unit.error_type and \
               'Z' not in err_unit.error_type:
                return ('X', err_unit.error_type, err_unit.error_source)
            if 'Z' in err_unit.error_type and \
               'X' not in err_unit.error_type and \
               'Y' not in err_unit.error_type:
                return ('Z', err_unit.error_type, err_unit.error_source)
        return ('MIXED', None, None)
    
    def _err_type_all(self) -> str:
        err_type_list = []
        for err_unit in self.err_unit_list:
            if 'X' in err_unit.error_type and \
               'Y' not in err_unit.error_type and \
               'Z' not in err_unit.error_type:
                err_type_list.append('X')
            if 'Z' in err_unit.error_type and \
               'X' not in err_unit.error_type and \
               'Y' not in err_unit.error_type:
                err_type_list.append('Z')
        if 'X' in err_type_list and \
           'Z' not in err_type_list:
            return 'X'
        elif 'X' in err_type_list and \
             'Z' in err_type_list:
            return 'BOTH'
        elif 'X' not in err_type_list and \
             'Z' in err_type_list:
            return 'Z'
        else:
            return 'MIXED'



class DEM():

    def __init__(self, stim_circuit: stim.Circuit, decompose: bool = False) -> None:
        self.dem_model = stim_circuit.detector_error_model(decompose_errors=decompose)
        self.num_dets = self.dem_model.num_detectors
        self.explained_err_list = stim_circuit\
                                .explain_detector_error_model_errors()
        self.fault_id_list, self.fault_unit_list = self._fault_units()
        self.fault_flipped_det_id_dict = self._fault_dict()
        self.x_faults, self.z_faults, \
            self.both_faults, self.mixed_faults = self._fault_cate()
        self.fault_id_dict: dict[int,Fault_Unit] = self._gen_id_dict()
        self.det_id_coords: dict[int,list[float]] = self.dem_model.get_detector_coordinates()
        self.post_selection_mask: list[int] = self.post_selected_det_ids()
        self.basis_det_assemble()
        pass

    def _gen_id_dict(self) -> dict[int,Fault_Unit]:
        ret_dict = {}
        for fault_unit in self.fault_unit_list:
            ret_dict[fault_unit.fault_id] = fault_unit
        return ret_dict

    def _fault_units(self) -> Tuple[list[int],list[Fault_Unit]]:
        fault_id_list = [i for i in range(len(self.explained_err_list))]
        fault_unit_list = []
        for explained_err, fault_id in zip(self.explained_err_list,
                                           fault_id_list):
            err_location_list = explained_err.circuit_error_locations
            dem_target_coord_list = explained_err.dem_error_terms
            dem_target_list = [dem_target_coord.dem_target for 
                               dem_target_coord in dem_target_coord_list]
            fault_unit_list.append(Fault_Unit(fault_id, 
                                              err_location_list,
                                              dem_target_list))
        return fault_id_list, fault_unit_list

    def _fault_dict(self) -> dict[frozenset,Tuple[bool,bool,float]]:
        ret_dict = {}
        for fault_unit in self.fault_unit_list:
            ret_dict[frozenset(fault_unit.flipped_detector_id_list)] = \
                (True, fault_unit.logical_flag, fault_unit.err_probability)
        return ret_dict

    def _fault_cate(self) -> Tuple[list[Fault_Unit],
                                   list[Fault_Unit],
                                   list[Fault_Unit],
                                   list[Fault_Unit]]:
        ret_x_faults = []
        ret_z_faults = []
        ret_both_faults = []
        ret_mixed_faults = []
        for fault_unit in self.fault_unit_list:
            if fault_unit.err_type_all == 'X':
                ret_x_faults.append(fault_unit)
            elif fault_unit.err_type_all == 'Z':
                ret_z_faults.append(fault_unit)
            elif fault_unit.err_type_all == 'BOTH':
                ret_x_faults.append(fault_unit)
                ret_z_faults.append(fault_unit)
                ret_both_faults.append(fault_unit)
            else:
                ret_mixed_faults.append(fault_unit)
        return ret_x_faults, ret_z_faults, ret_both_faults,\
               ret_mixed_faults
    
    def post_selected_det_ids(self) -> list[int]:
        ret_list = []
        for id in self.det_id_coords:
            if len(self.det_id_coords[id]) >= 4:
                ret_list.append(id)
        return ret_list
    
    def basis_det_assemble(self) -> None:
        self.x_det_ids = set()
        self.z_det_ids = set()
        for id in self.det_id_coords:
            if len(self.det_id_coords[id]) == 3:
                pos_x = self.det_id_coords[id][0]
                pos_y = self.det_id_coords[id][1]
                if (pos_x-pos_y) % 2 == 0:
                    self.x_det_ids.add(id)
                elif (pos_x-pos_y) % 2 == 1:
                    self.z_det_ids.add(id)
    
    def basis_cal(self, det_id) -> tuple[str, bool]:
        if det_id in self.x_det_ids:
            return ('X',True)
        elif det_id in self.z_det_ids:
            return ('Z',True)
        else:
            return ('', False)

    def prune_post_selected(self) -> stim.DetectorErrorModel:
        new_dem = stim.DetectorErrorModel()
        det_inst_list: list[stim.DemInstruction] = []
        err_inst_list: list[stim.DemInstruction] = []
        err_inst_basis: list[stim.DemInstruction] = []
        err_basis_dict: dict[frozenset[int],stim.DemInstruction] = {}
        err_inst_composite: list[stim.DemInstruction] = []
        err_inst_decomposed: list[stim.DemInstruction] = []
        for dem_inst in self.dem_model:
            if dem_inst.type == 'detector':
                det_id = dem_inst.targets_copy()[0].val
                det_coord = self.det_id_coords[det_id]
                det_coord_hint = []
                for pos in det_coord:
                    if pos % 1 == 0:
                        det_coord_hint.append(round(pos))
                    else:
                        det_coord_hint.append(pos)
                det_inst = stim.DemInstruction('detector',
                            det_coord_hint,
                            [stim.target_relative_detector_id(det_id)])
                det_inst_list.append(det_inst)
                # new_dem.append(det_inst)
            elif dem_inst.type == 'error':
                p_err = dem_inst.args_copy()[0]
                targets = dem_inst.targets_copy()
                is_post_selected = False
                for target in targets:
                    if target.is_relative_detector_id():
                        if target.val in self.post_selection_mask:
                            is_post_selected = True
                            break
                if is_post_selected == False:
                    err_inst = stim.DemInstruction('error',
                                    [p_err],targets)
                    err_inst_list.append(err_inst)
                    # new_dem.append(err_inst)
            else:
                raise ValueError('dem is wierd. this line is neither '
                                 'a detector or an error, consider '
                                 'flattening the dem.')

        while err_inst_list:
            err_inst = err_inst_list.pop()
            det_ids = []
            det_basis = []
            for target in err_inst.targets_copy():
                if target.is_relative_detector_id():
                    det_ids.append(target.val)
                    basis, _ = self.basis_cal(target.val)
                    if _ == False:
                        raise ValueError('det id not found')
                    det_basis.append(basis)
            if len(det_basis) == 1:
                err_inst_basis.append(err_inst)
                try:
                    err_basis_dict[frozenset(det_ids)]
                    raise ValueError('fault distance too small')
                except:
                    err_basis_dict[frozenset(det_ids)] = err_inst
            elif len(det_basis) == 2 and 'X' in det_basis \
                    and 'Z' not in det_basis:
                err_inst_basis.append(err_inst)
                try:
                    err_basis_dict[frozenset(det_ids)]
                    raise ValueError('fault distance too small')
                except:
                    err_basis_dict[frozenset(det_ids)] = err_inst
            elif len(det_basis) == 2 and 'Z' in det_basis \
                    and 'X' not in det_basis:
                err_inst_basis.append(err_inst)
                try:
                    err_basis_dict[frozenset(det_ids)]
                    raise ValueError('fault distance too small')
                except:
                    err_basis_dict[frozenset(det_ids)] = err_inst
            else:
                err_inst_composite.append(err_inst)
        self.err_basis_dict = err_basis_dict
        while err_inst_composite:
            err_inst = err_inst_composite.pop()
            det_ids = []
            for target in err_inst.targets_copy():
                if target.is_relative_detector_id():
                    det_ids.append(target.val)
            x_det_ids = [id for id in det_ids if id in self.x_det_ids]
            z_det_ids = [id for id in det_ids if id in self.z_det_ids]
            try:
                err_basis_1 = err_basis_dict[frozenset(x_det_ids)]
                err_basis_2 = err_basis_dict[frozenset(z_det_ids)]
                dem_target_list = err_basis_1.targets_copy() \
                                    + [stim.target_separator()] \
                                  + err_basis_2.targets_copy()
                err_inst_decomposed.append(
                                stim.DemInstruction('error',
                                            err_inst.args_copy(),
                                            dem_target_list))
            except:
                # print(x_det_ids,' x/z ',z_det_ids)
                # print('kind of wierd. Unknown error mechanism.')
                pass

        
        for err_inst in err_inst_basis:
            new_dem.append(err_inst)
        for err_inst in err_inst_decomposed:
            new_dem.append(err_inst)
        for det_inst in det_inst_list:
            new_dem.append(det_inst)
                
        self.det_id_connection: dict[int,list[frozenset[int]]] = {}
        for err_inst in err_inst_basis:
            det_ids = []
            for target in err_inst.targets_copy():
                if target.is_relative_detector_id():
                    det_ids.append(target.val)
            for det_id in det_ids:
                try:
                    self.det_id_connection[det_id].append(frozenset(det_ids))
                except:
                    self.det_id_connection[det_id] = [frozenset(det_ids)]

        # print(new_dem)
        return new_dem



def distance(coord_1: tuple[float], coord_2: tuple[float]) -> float:
    return np.max(np.abs(np.array(coord_1)-np.array(coord_2)))

def inv(coord: tuple[float]) -> tuple[float]:
    ret_list = [i for i in coord]
    ret_list[0] = -coord[0]
    ret_list[1] = -coord[1]
    return tuple(ret_list)