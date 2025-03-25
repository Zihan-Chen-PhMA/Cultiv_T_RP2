from __future__ import annotations

import numpy as np
from numpy import ndarray
from typing import Any, Union, TypeVar, Tuple
from warnings import warn
import stim


oddint = TypeVar('oddint')

class Operation():

    def __init__(self, operation: str, target: list[Qubit],
                 operation_id: int,
                 tick_id: int,
                 meas_id: Union[int, None],
                 clocktime: int,
                 arg: Union[float, None] = None
                 ) -> None:
        self.operation = operation
        self.target_id_list = [qubit.circuit_id for qubit in target]
        self.operation_id = operation_id
        self.tick_id = tick_id
        self.meas_id = meas_id
        self.clocktime = clocktime
        self.arg = None

class Measurement():

    def __init__(self, target: Qubit,
                 meas_id: Union[int, None],
                 clocktime: int
                 ) -> None:
        self.target = target
        self.meas_id = meas_id
        self.clocktime = clocktime

class Detector():

    def __init__(self, spacetime_coord: tuple[float], 
                 detector_id: int, 
                 list_meas_id: list[int] = None
                ) -> None:
        self.spacetime_coord = spacetime_coord
        self.list_meas_id = []
        self.detector_id = detector_id
        self.flag = False
        self.cz_flag = False
        pass

    def _append_meas_id(self, meas_id: int) -> None:
        self.list_meas_id.append(meas_id)
        pass

     

        
class Circuit_helper():
    # helper class for the stim package. 
    def __init__(self, qubit_collection: list[Qubit],
                 stim_filename: str
                 ) -> None:
        self.qubit_collection = qubit_collection
        self.stim_filename = stim_filename
        open(self.stim_filename, "w").close()

        self.time_temp: int = 0
        self.clock: list[int] = []
        self.time_meas_dict: dict[int,list[Measurement]] = {}
        self.tick_temp = 0
        self.tick_seq = [0]
        self.operation_id_temp = 0
        self.operation_id_seq = []
        self.operation_id_dict: dict[int,Operation] = {}
        self.meas_num = 0
        self.meas_id_list = []
        self.meas_oper_dict: dict[int,Operation] = {}
        self.detector_id_temp = 0
        self.detector_id_list = []
        self.detector_dict: dict[int, Detector] = {}
        self.spacetime_detector_dict: dict[Tuple,Detector] = {}

        self.id_qubit_dict: dict[int,Qubit] = {}
        self.num_qubits = None
        self.pos_qubit_dict: dict[Tuple,Qubit] = {}
        pass


    def initialize(self) -> None:

        self.num_qubits =  len(self.qubit_collection)
        circuit_id_list = [i for i in range(self.num_qubits)]
        for circuit_id, qubit in zip(circuit_id_list,
                                     self.qubit_collection):
            qubit.set_circuit_id(circuit_id)
            self.id_qubit_dict[circuit_id] = qubit
            self.pos_qubit_dict[qubit.pos] = qubit
        
        with open(self.stim_filename,'a') as f:
            for qubit in self.qubit_collection:
                f.write('QUBIT_COORDS'+qubit.pos_str+' '+str(qubit.circuit_id))
                f.write('\n')
        
        pass





    def reset(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None):
        if len(qubit_collection) == 0:
            return
        with open(self.stim_filename,'a') as f:
            f.write('R')
            for qubit in qubit_collection:
                oper = Operation('R', [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                self.operation_id_seq.append(self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')

    def reset_X(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None):
        if len(qubit_collection) == 0:
            return
        with open(self.stim_filename,'a') as f:
            f.write('RX')
            for qubit in qubit_collection:
                oper = Operation('RX', [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                self.operation_id_seq.append(self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'Z_ERROR':
                warn('Not an Z_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('Z_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')

    def reset_Y(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None):
        with open(self.stim_filename,'a') as f:
            f.write('RY')
            for qubit in qubit_collection:
                oper = Operation('RY', [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                self.operation_id_seq.append(self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'DEPOLARIZE1':
                warn('Not an DEPOLARIZE1 ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('DEPOLARIZE1')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
    
    def MPP(self, meas_supps: list[tuple[Qubit]], 
            meas_opers: list[tuple[str]], 
            qubit_register: list[Qubit]) -> None:
        meas_oper_list = []
        with open(self.stim_filename,'a') as f:
            f.write('MPP')
            for meas_supp, meas_oper, register in zip(meas_supps,
                                    meas_opers, qubit_register):
                for i in range(len(meas_supp)):
                    if i != 0:
                        f.write('*')
                    else:
                        f.write(' ')
                    f.write(meas_oper[i] + 
                            str(meas_supp[i].circuit_id))
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                register.meas_id_seq.append(meas_id)
                register.time_meas_id_seq.append((self.time_temp, meas_id))
                register.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('M', [register], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                register.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
                
                




    def single_qubit_gate(self, gate_name: str, qubit_collection: list[Qubit],
                          noise: Union[list[str,float],None] = None):
        if qubit_collection == None:
            return None
        if len(qubit_collection) == 0:
            return None
        with open(self.stim_filename,'a') as f:
            f.write(gate_name)
            for qubit in qubit_collection:
                oper = Operation(gate_name, [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    if qubit != None:
                        f.write(' ')
                        f.write(str(qubit.circuit_id))
                f.write('\n')
    
    def single_qubit_noise(self, qubit_collection: list[Qubit],
                           noise: Union[list[str,float],None] = None):
        if qubit_collection == None:
            return None
        if len(qubit_collection) == 0:
            return None
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    if qubit != None:
                        f.write(' ')
                        f.write(str(qubit.circuit_id))
                f.write('\n')
        
    
    def two_qubit_gate(self, gate_name: str, 
                       control_qubit_collection: list[Qubit],
                       target_qubit_collection: list[Qubit],
                       noise: Union[list[str,float],None] = None):
        if control_qubit_collection == None:
            return None
        if len(control_qubit_collection) == 0:
            return None
        if len(control_qubit_collection) != len(target_qubit_collection):
            raise ValueError('control-target mismatch')
        with open(self.stim_filename,'a') as f:
            f.write(gate_name)
            for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                              target_qubit_collection):
                if ctrl_qubit != None and targ_qubit != None:
                    oper = Operation(gate_name, [ctrl_qubit, targ_qubit], 
                                     self.operation_id_temp,
                                     self.tick_temp, None, self.time_temp)
                    self.operation_id_dict[self.operation_id_temp] = oper
                    self.operation_id_seq.append(self.operation_id_temp)
                    ctrl_qubit.append_operation(self.tick_temp, 
                                                self.operation_id_temp)
                    targ_qubit.append_operation(self.tick_temp,
                                                self.operation_id_temp)
                    f.write(' ')
                    f.write(str(ctrl_qubit.circuit_id))
                    f.write(' ')
                    f.write(str(targ_qubit.circuit_id))
                    self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                                  target_qubit_collection):
                    if ctrl_qubit != None and targ_qubit != None:         
                        f.write(' ')
                        f.write(str(ctrl_qubit.circuit_id))
                        f.write(' ')
                        f.write(str(targ_qubit.circuit_id))
                f.write('\n')

        pass

    def two_qubit_noise(self, 
                       control_qubit_collection: list[Qubit],
                       target_qubit_collection: list[Qubit],
                       noise: Union[list[str,float],None] = None):
        if control_qubit_collection == None:
            return None
        if len(control_qubit_collection) == 0:
            return None
        if len(control_qubit_collection) != len(target_qubit_collection):
            raise ValueError('control-target mismatch')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                                  target_qubit_collection):
                    if ctrl_qubit != None and targ_qubit != None:         
                        f.write(' ')
                        f.write(str(ctrl_qubit.circuit_id))
                        f.write(' ')
                        f.write(str(targ_qubit.circuit_id))
                f.write('\n')

        pass

    def measure_reset(self, qubit_collection: list[Qubit],
                      noise: Union[list[str,float],None] = None) -> None:
        if len(qubit_collection) == 0:
            return
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('MR')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_meas_id_seq.append((self.time_temp, meas_id))
                qubit.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('MR', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        # self.time_temp = self.time_temp + 1
        pass

    def flag_measure(self, qubit_collection: list[Qubit],
                     stage: int,
                     noise: Union[list[str,float],None] = None) -> None:
        if len(qubit_collection) == 0:
            return
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('M')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_flag_meas_id_seq.append((stage, meas_id))
                qubit.time_flag_meas_id_dict[stage] = meas_id
                meas_oper = Operation('M', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list

    def flag_measure_X(self, qubit_collection: list[Qubit], 
                       stage: int,
                       noise: Union[list[str,float],None] = None) -> None:
        if len(qubit_collection) == 0:
            return
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'Z_ERROR':
                warn('Not an Z_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('Z_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('MX')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_flag_meas_id_seq.append((stage, meas_id))
                qubit.time_flag_meas_id_dict[stage] = meas_id
                meas_oper = Operation('MX', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        pass
        

    def measure(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None,
                not_qubits: list[Qubit] = []) -> None:
        if len(qubit_collection) == 0:
            return
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('M')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_meas_id_seq.append((self.time_temp, meas_id))
                qubit.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('M', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                if qubit in not_qubits:
                    f.write('!')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        pass

    def measure_X(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None,
                not_qubits: list[Qubit] = []) -> None:
        if len(qubit_collection) == 0:
            return
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'Z_ERROR':
                warn('Not an Z_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('Z_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('MX')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_meas_id_seq.append((self.time_temp, meas_id))
                qubit.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('MX', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                if qubit in not_qubits:
                    f.write('!')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        pass

    def clock_plus_one(self) -> None:
        self.time_temp = self.time_temp + 1
    
    def feedforward_X(self, meas_qubit_list: list[Qubit],
                      X_correction_qubit_list: list[Qubit],
                      time_list: list[int]) -> None:
        if len(meas_qubit_list) == 0:
            raise ValueError('measurements unspecified')
        if len(meas_qubit_list) != len(X_correction_qubit_list) or \
            len(meas_qubit_list) != len(time_list):
            raise ValueError('input array size mismatch')
        with open(self.stim_filename,'a') as f:
            f.write('CX') 
            for meas_qubit, corr_qubit, time in zip(
                                meas_qubit_list,
                                X_correction_qubit_list,
                                time_list):
                # search for meas_id at time_front:
                index = -1
                for i, time_meas_id_tuple in enumerate(meas_qubit.time_meas_id_seq):
                    if time == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('undefined time step')
                meas_id = meas_qubit.time_meas_id_seq[index][-1]
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(' rec')
                f.write(str([neg_pos]))
                f.write(' ')
                f.write(str(corr_qubit.circuit_id))
            f.write('\n')
    
    def feedforward_pauli(self, pauli: str, meas_qubit_list: list[Qubit],
                 corr_qubit_list: list[Qubit],
                 time_or_stage_list: list[tuple[int,str]]) -> None:
        if len(meas_qubit_list) == 0:
            raise ValueError('measurements unspecified')
        if len(meas_qubit_list) != len(corr_qubit_list) or \
            len(meas_qubit_list) != len(time_or_stage_list):
            raise ValueError('input array size mismatch')
        with open(self.stim_filename,'a') as f:
            f.write('C')
            f.write(pauli) 
            for meas_qubit, corr_qubit, time_or_stage in zip(
                                meas_qubit_list,
                                corr_qubit_list,
                                time_or_stage_list):
                # search for meas_id at time_front:
                if time_or_stage[-1] == 'time':
                    time = time_or_stage[0]
                    try:
                        meas_id = meas_qubit.time_meas_id_dict[time]
                    except:
                        raise ValueError('undefined time step')
                elif time_or_stage[-1] == 'stage':
                    stage = time_or_stage[0]
                    try:
                        meas_id = meas_qubit.time_flag_meas_id_dict[stage]
                    except:
                        raise ValueError('undefined flag stage')
                
                assert meas_id
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(' rec')
                f.write(str([neg_pos]))
                f.write(' ')
                f.write(str(corr_qubit.circuit_id))
            f.write('\n')
        

    def detector(self, qubit_front_list: list[Qubit], time_front: int,
                 qubit_back_list: list[Qubit], time_back: int, 
                 spacetime_coord: tuple[float],
                 flag: bool = False,
                 cz_flag: bool = False) -> None:
        if len(qubit_front_list) + len(qubit_back_list) == 0:
            raise ValueError('operation not supported on anything')
        ret_detector = Detector(spacetime_coord, 
                                self.detector_id_temp)
        if flag == True:
            ret_detector.flag = True
        if cz_flag == True:
            ret_detector.cz_flag = True
        self.detector_id_list.append(self.detector_id_temp)
        spacetime_coord_str = pos_str_conversion_general([coord for
                                                    coord in spacetime_coord])
        with open(self.stim_filename,'a') as f:
            f.write('DETECTOR')
            f.write(spacetime_coord_str)
            for qubit in qubit_front_list:
                if qubit == None:
                    continue
                f.write(' rec')
                try:
                    meas_id = qubit.time_meas_id_dict[time_front]
                except:
                    raise ValueError('undefined time step')
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
            for qubit in qubit_back_list:
                if qubit == None:
                    continue
                f.write(' rec')
                try:
                    meas_id = qubit.time_meas_id_dict[time_back]
                except:
                    raise ValueError('undefined time step')
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
            f.write('\n')
        self.detector_dict[self.detector_id_temp] = ret_detector
        self.spacetime_detector_dict[ret_detector.spacetime_coord] = \
                                                    ret_detector
        self.detector_id_temp = self.detector_id_temp + 1
        pass
    
    def flag_detector(self, qubit: Qubit,
                      stage: int, spacetime_coord) -> None:
        ret_detector = Detector(spacetime_coord, 
                                self.detector_id_temp)
        
        self.detector_id_list.append(self.detector_id_temp)
        spacetime_coord_str = pos_str_conversion_general([coord for
                                                    coord in spacetime_coord])
        with open(self.stim_filename,'a') as f:
            f.write('DETECTOR')
            f.write(spacetime_coord_str)
            f.write(' rec')
            try:
                meas_id = qubit.time_flag_meas_id_dict[stage]
            except:
                raise ValueError('undefined time step')
            ret_detector._append_meas_id(meas_id)
            neg_pos = (self.meas_id_list.index(meas_id) 
                        -len(self.meas_id_list))
            f.write(str([neg_pos]))

            f.write('\n')
        self.detector_dict[self.detector_id_temp] = ret_detector
        self.spacetime_detector_dict[ret_detector.spacetime_coord] = \
                                                    ret_detector
        self.detector_id_temp = self.detector_id_temp + 1

    def detector_qt_list(self, 
                 qubit_list: list[Qubit], 
                 time_list: list[int],
                 spacetime_coord: tuple[float],
                 flag: bool = False,
                 cz_flag: bool = False) -> None:
        if len(qubit_list) == 0:
            raise ValueError('detector on null set')
        if len(qubit_list) != len(time_list):
            raise ValueError('length mismatch')
        ret_detector = Detector(spacetime_coord, 
                                self.detector_id_temp)
        if flag == True:
            ret_detector.flag = True
        if cz_flag == True:
            ret_detector.cz_flag = True
        self.detector_id_list.append(self.detector_id_temp)
        spacetime_coord_str = pos_str_conversion_general([coord for
                                                    coord in spacetime_coord])
        with open(self.stim_filename,'a') as f:
            f.write('DETECTOR')
            f.write(spacetime_coord_str)
            for qubit, time_q in zip(qubit_list, time_list):
                if qubit == None:
                    continue
                f.write(' rec')
                try:
                    meas_id = qubit.time_meas_id_dict[time_q]
                except:
                    raise ValueError('undefined time step')
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))

            f.write('\n')
        self.detector_dict[self.detector_id_temp] = ret_detector
        self.spacetime_detector_dict[ret_detector.spacetime_coord] = \
                                                    ret_detector
        self.detector_id_temp = self.detector_id_temp + 1
        pass


    def observable(self, qubit_list: list[Qubit], time: int) -> None:
        with open(self.stim_filename,'a') as f:
            f.write('OBSERVABLE_INCLUDE')
            f.write('(0)')
            for qubit in qubit_list:
                f.write(' rec')
                # search for meas_id at time_front:
                index = -1
                for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                    if time == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('undefined time step')
                meas_id = qubit.time_meas_id_seq[index][-1]
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
            f.write('\n')
        pass


    def tick(self) -> None:
        with open(self.stim_filename,'a') as f:
            f.write('TICK\n')
        self.tick_temp = self.tick_temp + 1
        pass
    
    
    

       

class Qubit():
    
    def __init__(self, id: int, pos: tuple[float, float], tribe: str = None,
                 tag: str = None
                 ) -> None:
        self.id = id
        self.pos = pos
        self.pos_str = pos_str_conversion(self.pos)
        self.tribe = tribe     # 'D', 'X', 'Z', 'FLAG'
        self.neighbor: dict[str,Union[Qubit,None]] = {}
        self.check = {}
        self.css = None
        self.tag = tag
        # self.oper_seq = []
        self.meas_id_seq = []
        self.time_meas_id_seq = []
        self.time_meas_id_dict: dict[int,int] = {}
        self.time_flag_meas_id_seq = []
        self.time_flag_meas_id_dict: dict[int,int] = {}
        self.circuit_id = None
        self.tick_oper_id_dict: dict[int, int] = {}
    
    def append_neighbor(self, qubit: Union[Qubit,None], rel_pos_str: str):
        self.neighbor[rel_pos_str] = qubit

    def append_check(self, qubit: Qubit, single_check: str):
        self.check[pos_str_conversion(qubit.pos)] = [qubit, single_check]
        pass

    def set_circuit_id(self, circuit_id: int) -> None:
        self.circuit_id = circuit_id

    def append_operation(self, tick: int, oper_id: int) -> None:
        self.tick_oper_id_dict[tick] = oper_id
        pass





def pos_str_conversion(pos: tuple[Union[int,oddint], 
                                  Union[int,oddint]]) -> str:
    pos_x = pos[0]
    pos_y = pos[1]
    if abs(pos_x-round(pos_x)) < 10**(-8):
        pos_x_str = "{:.0f}".format(pos_x)
    elif abs(abs(pos_x-round(pos_x))-0.5) < 10**(-8):
        pos_x_str = "{:.1f}".format(pos_x)
    else:
        return ValueError('pos_x has a strange value')
    if abs(pos_y-round(pos_y)) < 10**(-8):
        pos_y_str = "{:.0f}".format(pos_y)
    elif abs(abs(pos_y-round(pos_y))-0.5) < 10**(-8):
        pos_y_str = "{:.1f}".format(pos_y)
    else:
        return ValueError('pos_y has a strange value')
    pos_str = '(' + pos_x_str + ',' + pos_y_str + ')'
    return pos_str


def pos_str_conversion_general(pos: list[Union[int, oddint]]) -> str:
    ret_str = '('
    for coord, index in zip(pos,range(len(pos))):
        if abs(coord-round(coord)) < 10**(-8):
            coord_str = "{:.0f}".format(coord)
        elif abs(abs(coord-round(coord))-0.5) < 10**(-8):
            coord_str = "{:.1f}".format(coord)
        else:
            return ValueError('coord has a strange value')
        ret_str = ret_str + coord_str
        if index != len(pos)-1:
            ret_str = ret_str + ','
    ret_str = ret_str + ')'
    return ret_str



def grid_network(lx: float, ly: float,
                 x_gap: float, y_gap: float,
                 x_offset: float, y_offset: float, 
                 nodes: list[Qubit]) -> Tuple[dict[str, Qubit],
                                           dict[str, tuple[float, float]]]: 
    if len(nodes) != lx*ly:
        raise ValueError('length of the list of the nodes is incompatible')
    network = {}
    network_coord = {}
    for i in range(lx):
        for j in range(ly):
            pos_x = i*x_gap + x_offset
            pos_y = j*y_gap + y_offset
            pos_str = pos_str_conversion((pos_x, pos_y))
            network[pos_str] = nodes[i*ly + j]
            network_coord[pos_str] = (pos_x,pos_y)
    return [network, network_coord]








class RP_surface_code_fresh():


    def __init__(self, distance: oddint, code_name: Union[str,None] = None) -> None:
        if (distance%2) != 1:
            raise ValueError('distance is not an odd number') 
        self.distance = distance
        self.x_displacement = self.distance // 2    # make sure centre qubit is at (0,0)
        self.y_displacement = self.distance // 2
        self.code_name = code_name
        self.num_qubits = distance**2
        self.id_collection = [i for i in range(self.num_qubits)]
        self.data_qubits_collection = self.initialize_data_qubits()
        self.data_qubits_network: dict[tuple,Qubit] = {}
        self.generate_data_qubit_network()
        self.z_check_collection = self.initialize_z_checks()
        self.generate_z_check_network()
        self.num_z_checks = len(self.z_check_collection)
        self.x_check_collection = self.initilize_x_checks()
        self.num_x_checks = len(self.x_check_collection)
        self.generate_x_check_network()
        self.qubit_network: dict[tuple,Qubit] = {}
        self.logic_z_collection = self._logic_z_assemble()
        self.logic_x_collection = self._logic_x_assemble()
        self.logic_y_supp, self.logic_y_oper = self._logic_y_assemble()
        self._qubit_displacement()
        for qubit in self.data_qubits_collection \
                     + self.z_check_collection \
                     + self.x_check_collection:
            self.qubit_network[qubit.pos] = qubit
        pass


    def _logic_z_assemble(self) -> list[Qubit]:
        x_pos = (self.distance-1)/2
        ret_z_list = []
        for qubit in self.data_qubits_collection:
            if x_pos == qubit.pos[0]:
                ret_z_list.append(qubit)
        return ret_z_list

    def _logic_x_assemble(self) -> list[Qubit]:
        y_pos = (self.distance-1)/2
        ret_x_list = []
        for qubit in self.data_qubits_collection:
            if y_pos == qubit.pos[1]:
                ret_x_list.append(qubit)
        return ret_x_list
    
    def _logic_y_assemble(self) -> tuple[list[Qubit],
                                         list[str]]:
        supps = []
        opers = []
        for qubit in self.logic_z_collection:
            supps.append(qubit)
            if qubit in self.logic_x_collection:
                opers.append('Y')
            else:
                opers.append('Z')
        for qubit in self.logic_x_collection:
            if qubit not in self.logic_z_collection:
                supps.append(qubit)
                opers.append('X')
        return (supps, opers)
            


    def initialize_data_qubits(self) -> list[Qubit]:
        data_qubits_collection = []
        for i in range(self.distance):
            for j in range(self.distance):
                order = i*self.distance + j
                pos = (i,j)
                id = self.id_collection[order]
                data_qubits_collection.append(Qubit(id, pos, 'D'))
        return data_qubits_collection
    
    def _qubit_displacement(self) -> None:
        for qubit in self.data_qubits_collection \
                     + self.x_check_collection \
                     + self.z_check_collection:
            qubit.pos = (qubit.pos[0]-self.x_displacement,
                         qubit.pos[1]-self.y_displacement)

    def initilize_x_checks(self) -> list[Qubit]: # Unfinished. 
        x_check_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_left_bd_collection = []
        id_upper_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 0.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_A')
                x_check_collection.append(x_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 1.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_B')
                x_check_collection.append(x_check_temp)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_left_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = -0.5
            pos_y = 2*j + 1.5
            x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'left_bd')
            x_check_collection.append(x_check_temp)
        # print(id_collection)
        for j in range((self.distance//2) + 1):
            id_index = id_index + 1
            id_upper_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = 2*j - 0.5
            pos_y = -0.5
            if j > 0:
                x_check_temp = Qubit(id_index,(pos_x, pos_y), 'X', 'upper_bd')
            else:
                x_check_temp = Qubit(id_index,(pos_x,pos_y),'X','corner')
            x_check_collection.append(x_check_temp)
        
        return x_check_collection

    
    def initialize_z_checks(self) -> list[Qubit]:
        z_check_collection = []
        id_upper_bd_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_left_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range((self.distance//2) + 1):
            id_index = id_index + 1
            id_upper_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = i*2 + 0.5
            pos_y = -0.5
            if i < (self.distance//2):
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'upper_bd')
            else:
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'corner')
            
            z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 1.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_A')
                z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 0.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_B')
                z_check_collection.append(z_check_temp)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_left_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = -0.5
            pos_y = 2*j + 0.5
            z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'left_bd')
            z_check_collection.append(z_check_temp)
        return z_check_collection
    



    def generate_data_qubit_network(self):

        self.data_qubits_network: dict[tuple,Qubit] = {}
        for data in self.data_qubits_collection:
            self.data_qubits_network[data.pos] = data
        pass


    def generate_x_check_network(self):
        for check in self.x_check_collection:
            if check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                if check.pos[1] == self.distance-1.5 \
                    and check.tag=='sub_B':
                    check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                    check.append_neighbor(self.data_qubits_network[se_pos],'5')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'4')
                elif check.pos[0] == self.distance-2.5 \
                    and check.tag == 'sub_A':
                    check.append_neighbor(self.data_qubits_network[sw_pos],'2')
                    check.append_neighbor(self.data_qubits_network[se_pos],'5')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'3')
                else:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'2')
                    check.append_neighbor(self.data_qubits_network[se_pos],'4')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'3')
                check.append_check(self.data_qubits_network[sw_pos],'X')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[nw_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
            elif check.tag == 'left_bd':
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                sw_pos = (self.distance-1, self.distance-check.pos[1]-1.5)
                nw_pos = (self.distance-1, self.distance-check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[sw_pos],'1')
                check.append_neighbor(self.data_qubits_network[se_pos],'4')
                check.append_neighbor(self.data_qubits_network[nw_pos],'2')
                check.append_neighbor(self.data_qubits_network[ne_pos],'5')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
                check.append_check(self.data_qubits_network[sw_pos],'X')
                check.append_check(self.data_qubits_network[nw_pos],'X')
            elif check.tag == 'upper_bd':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (self.distance-check.pos[0]-0.5, self.distance-1)
                ne_pos = (self.distance-check.pos[0]-1.5, self.distance-1)
                check.append_neighbor(self.data_qubits_network[sw_pos],'5')
                check.append_neighbor(self.data_qubits_network[se_pos],'3')
                check.append_neighbor(self.data_qubits_network[nw_pos],'2')
                check.append_neighbor(self.data_qubits_network[ne_pos],'1')
                check.append_check(self.data_qubits_network[sw_pos],'X')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[nw_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
            elif check.tag == 'corner':
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                ne_pos = (self.distance-1,self.distance-1)
                check.append_neighbor(self.data_qubits_network[se_pos],'3')
                check.append_neighbor(self.data_qubits_network[ne_pos],'1')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
                pass
        pass


    def generate_z_check_network(self):
        for check in self.z_check_collection:
            if check.tag == 'upper_bd':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (self.distance-check.pos[0]-0.5, self.distance-1)
                ne_pos = (self.distance-check.pos[0]-1.5, self.distance-1)
                if check.pos[0] == 0.5:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'5')
                    check.append_neighbor(self.data_qubits_network[se_pos],'4')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'3')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'1')
                else:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'5')
                    check.append_neighbor(self.data_qubits_network[se_pos],'4')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'2')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'1')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[se_pos],'Z')
                check.append_check(self.data_qubits_network[nw_pos],'Z')
                check.append_check(self.data_qubits_network[ne_pos],'Z')
            elif check.tag == 'left_bd':
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                sw_pos = (self.distance-1, self.distance-check.pos[1]-1.5)
                nw_pos = (self.distance-1, self.distance-check.pos[1]-0.5)
                if check.pos[1] == 0.5:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[se_pos],'3')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'2')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'4')
                else:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[se_pos],'3')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'2')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'5')
                check.append_check(self.data_qubits_network[se_pos],'Z')
                check.append_check(self.data_qubits_network[ne_pos],'Z')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[nw_pos],'Z')
            elif check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                if check.tag == 'sub_B' \
                    and check.pos[0] == self.distance-1.5:
                    if check.pos[1] == 0.5:
                        check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                        check.append_neighbor(self.data_qubits_network[se_pos],'5')
                        check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                        check.append_neighbor(self.data_qubits_network[ne_pos],'2')
                    else:
                        check.append_neighbor(self.data_qubits_network[sw_pos],'4')
                        check.append_neighbor(self.data_qubits_network[se_pos],'5')
                        check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                        check.append_neighbor(self.data_qubits_network[ne_pos],'3')
                elif check.tag == 'sub_B' \
                    and check.pos[1] == self.distance-2.5:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                    check.append_neighbor(self.data_qubits_network[se_pos],'5')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'2')
                else:
                    check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                    check.append_neighbor(self.data_qubits_network[se_pos],'4')
                    check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                    check.append_neighbor(self.data_qubits_network[ne_pos],'2')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[se_pos],'Z')
                check.append_check(self.data_qubits_network[nw_pos],'Z')
                check.append_check(self.data_qubits_network[ne_pos],'Z')
            elif check.tag == 'corner':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (0,self.distance-1)
                check.append_neighbor(self.data_qubits_network[sw_pos],'4')
                check.append_neighbor(self.data_qubits_network[se_pos],'5')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[se_pos],'Z')
        pass



class Rotated_surface_code():


    def __init__(self, distance: oddint, code_name
                    :Union[str,None] = None) -> None:
        if (distance%2) != 1:
            raise ValueError('distance is not an odd number') 
        self.distance = distance
        self.x_displacement = self.distance // 2
        self.y_displacement = self.distance // 2
        self.code_name = code_name
        self.num_qubits = distance**2
        self.id_collection = [i for i in range(self.num_qubits)]
        self.data_qubits_collection = self.initialize_data_qubits()
        self.data_qubits_network: dict[tuple,Qubit] = {}
        self.generate_data_qubit_network()
        self.z_check_collection = self.initialize_z_checks()
        self.generate_z_check_network()
        self.num_z_checks = len(self.z_check_collection)
        self.x_check_collection = self.initilize_x_checks()
        self.num_x_checks = len(self.x_check_collection)
        self.generate_x_check_network()
        self.qubit_network: dict[tuple,Qubit] = {}
        for qubit in self.data_qubits_collection \
                     + self.z_check_collection \
                     + self.x_check_collection:
            self.qubit_network[qubit.pos] = qubit
        self.logic_z_collection = self._logic_z_assemble()
        self.logic_x_collection = self._logic_x_assemble()
        self.logic_y_supp, self.logic_y_oper = self._logic_y_assemble()
        self._qubit_displacement()
        for qubit in self.data_qubits_collection \
                     + self.z_check_collection \
                     + self.x_check_collection:
            self.qubit_network[qubit.pos] = qubit
        pass


    def _logic_z_assemble(self) -> list[Qubit]:
        x_pos = (self.distance-1)/2
        ret_z_list = []
        for qubit in self.data_qubits_collection:
            if x_pos == qubit.pos[0]:
                ret_z_list.append(qubit)
        return ret_z_list

    def _logic_x_assemble(self) -> list[Qubit]:
        y_pos = (self.distance-1)/2
        ret_x_list = []
        for qubit in self.data_qubits_collection:
            if y_pos == qubit.pos[1]:
                ret_x_list.append(qubit)
        return ret_x_list
    
    def _logic_y_assemble(self) -> tuple[list[Qubit],
                                         list[str]]:
        supps = []
        opers = []
        for qubit in self.logic_z_collection:
            supps.append(qubit)
            if qubit in self.logic_x_collection:
                opers.append('Y')
            else:
                opers.append('Z')
        for qubit in self.logic_x_collection:
            if qubit not in self.logic_z_collection:
                supps.append(qubit)
                opers.append('X')
        return (supps, opers)
    
    def _qubit_displacement(self) -> None:
        for qubit in self.data_qubits_collection \
                     + self.x_check_collection \
                     + self.z_check_collection:
            qubit.pos = (qubit.pos[0]-self.x_displacement,
                         qubit.pos[1]-self.y_displacement)
    
    def qubit_displacement(self, x_displacement,
                           y_displacement) -> None:
        for qubit in self.data_qubits_collection \
                     + self.x_check_collection \
                     + self.z_check_collection:
            qubit.pos = (qubit.pos[0]-x_displacement,
                         qubit.pos[1]-y_displacement)
    
    def logic_x_selection(self, y_pos) -> None:
        self.logic_x_collection = []
        for qubit in self.data_qubits_collection:
            if y_pos == qubit.pos[1]:
                self.logic_x_collection.append(qubit)
        return self.logic_x_collection
    
    def logic_z_selection(self, x_pos) -> list[Qubit]:
        self.logic_z_collection = []
        for qubit in self.data_qubits_collection:
            if x_pos == qubit.pos[0]:
                self.logic_z_collection.append(qubit)
        return self.logic_z_collection
 

    def initialize_data_qubits(self) -> list[Qubit]:
        data_qubits_collection = []
        for i in range(self.distance):
            for j in range(self.distance):
                order = i*self.distance + j
                pos = (i,j)
                id = self.id_collection[order]
                data_qubits_collection.append(Qubit(id, pos, 'D'))
        return data_qubits_collection
    
    def initilize_x_checks(self) -> list[Qubit]: # Unfinished. 
        x_check_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_left_bd_collection = []
        id_right_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 0.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_A')
                x_check_collection.append(x_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 1.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_B')
                x_check_collection.append(x_check_temp)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_left_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = -0.5
            pos_y = 2*j + 1.5
            x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'left_bd')
            x_check_collection.append(x_check_temp)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_right_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = self.distance - 0.5
            pos_y = 2*j + 0.5
            x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'right_bd')
            x_check_collection.append(x_check_temp)
        
        return x_check_collection

    
    def initialize_z_checks(self) -> list[Qubit]:
        z_check_collection = []
        id_upper_bd_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_lower_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range((self.distance//2)):
            id_index = id_index + 1
            id_upper_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = i*2 + 0.5
            pos_y = -0.5
            z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'upper_bd')
            z_check_collection.append(z_check_temp)
        for i in range((self.distance//2)):
            id_index = id_index + 1
            id_lower_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = i*2 + 1.5
            pos_y = self.distance - 0.5
            z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'lower_bd')
            z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 1.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_A')
                z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 0.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_B')
                z_check_collection.append(z_check_temp)
        return z_check_collection
    



    def generate_data_qubit_network(self):

        self.data_qubits_network: dict[tuple,Qubit] = {}
        for data in self.data_qubits_collection:
            self.data_qubits_network[data.pos] = data
        pass


    def generate_x_check_network(self):
        for check in self.x_check_collection:
            if check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[sw_pos],'2')
                check.append_neighbor(self.data_qubits_network[se_pos],'4')
                check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                check.append_neighbor(self.data_qubits_network[ne_pos],'3')
                check.append_check(self.data_qubits_network[sw_pos],'X')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[nw_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
            elif check.tag == 'left_bd':
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[se_pos],'4')
                check.append_neighbor(self.data_qubits_network[ne_pos],'3')
                check.append_check(self.data_qubits_network[se_pos],'X')
                check.append_check(self.data_qubits_network[ne_pos],'X')
            elif check.tag == 'right_bd':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[sw_pos],'2')
                check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                check.append_check(self.data_qubits_network[sw_pos],'X')
                check.append_check(self.data_qubits_network[nw_pos],'X')
                pass
        pass


    def generate_z_check_network(self):
        for check in self.z_check_collection:
            if check.tag == 'upper_bd':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                check.append_neighbor(self.data_qubits_network[se_pos],'4')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[se_pos],'Z')
            elif check.tag == 'lower_bd':
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                check.append_neighbor(self.data_qubits_network[ne_pos],'2')
                check.append_check(self.data_qubits_network[ne_pos],'Z')
                check.append_check(self.data_qubits_network[nw_pos],'Z')
            elif check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                check.append_neighbor(self.data_qubits_network[sw_pos],'3')
                check.append_neighbor(self.data_qubits_network[se_pos],'4')
                check.append_neighbor(self.data_qubits_network[nw_pos],'1')
                check.append_neighbor(self.data_qubits_network[ne_pos],'2')
                check.append_check(self.data_qubits_network[sw_pos],'Z')
                check.append_check(self.data_qubits_network[se_pos],'Z')
                check.append_check(self.data_qubits_network[nw_pos],'Z')
                check.append_check(self.data_qubits_network[ne_pos],'Z')

        pass



