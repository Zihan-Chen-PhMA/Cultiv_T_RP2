from __future__ import annotations
from RP_surface_codes import *
from typing import Callable
import copy
import pymatching
      


class Meta_Circuit():

    def __init__(self, code_list: list[Union[RP_surface_code_fresh,
                                             Rotated_surface_code]], 
                 circuit_filename: str,
                 noise_probability: float = 0, 
                 control_functional: Callable = None) -> None:
        self.code_list = code_list
        self.code_name_dict: dict[str,Union[RP_surface_code_fresh,
                                            Rotated_surface_code]] = {}
        for code in self.code_list:
            self.code_name_dict[code.code_name] = code
        self.qubit_collection: list[Qubit] = []
        self.qubit_network: dict[tuple,Qubit] = {}
        self._collect_qubits()
        self.code_name_SE_instruction: dict[str,SE_instruction_set] = {}
        for code_name in self.code_name_dict:
            self.code_name_SE_instruction[code_name] = self._SE_compile(code_name)

        self.circuit_filename = circuit_filename
        self.circ = Circuit_helper(self.qubit_collection,
                            circuit_filename)
        self.noise_probability = noise_probability
        self.growth_time = -1
        self.circ.initialize()
        self.logic_qubits = []
        self.active_qubit_ids = []
        self.used_qubit_ids = []
        self.used_qubits: list[Qubit] = []
        self.total_qubit_circ_ids = []
        self.circ_id_qubit_dict: dict[int,Qubit] = {}
        self.stage = 0
        self.stage_list = []
        self._circ_id_init()

    def _idle_depo(self, cleaness: bool = False) -> None:
        idle_qubits = []
        if cleaness == True:
            self.used_qubits = []
            self.used_qubit_ids = []
            return
        for qubit in self.used_qubits:
            self.used_qubit_ids.append(qubit.circuit_id)
        for i in self.active_qubit_ids:
            if i not in self.used_qubit_ids:
                idle_qubits.append(self.circ_id_qubit_dict[i])
        self.used_qubits = []
        self.used_qubit_ids = []
        self.circ.single_qubit_noise(idle_qubits,
                                     noise=['DEPOLARIZE1',self.noise_probability])


    def _add_to_active(self,reset_qubit_list:list[Qubit]) -> None:
        for qubit in reset_qubit_list:
            if qubit.circuit_id not in self.active_qubit_ids:
                self.active_qubit_ids.append(qubit.circuit_id)
            else:
                raise ValueError('Qubit already active.' 
                                 'You may have reset an active qubit.')
        pass

    def _pop_from_active(self,meas_qubit_list:list[Qubit]) -> None:
        for qubit in meas_qubit_list:
            try:
                self.active_qubit_ids.remove(qubit.circuit_id)
            except:
                raise ValueError('Measured an inactive qubit.')            
    

    def _circ_id_init(self) -> None:
        for qubit in self.qubit_collection:
            self.total_qubit_circ_ids.append(qubit.circuit_id)
            self.circ_id_qubit_dict[qubit.circuit_id] = qubit

    def _collect_qubits(self) -> None:
        self.qubit_collection = []
        qubit_y_pos = set()
        qubit_x_pos = set()
        self.flag_network = {}
        id_temp = 0
        for code in self.code_list:
            for qubit in code.qubit_network.values():
                qubit_x_pos.add(qubit.pos[0])
                qubit_y_pos.add(qubit.pos[1])
                if qubit.pos not in self.qubit_network:
                    qubit_new = Qubit(id_temp,qubit.pos)
                    self.qubit_collection.append(qubit_new)
                    self.qubit_network[qubit.pos] = qubit_new
        
        max_y_abs = round(max(np.abs(list(qubit_y_pos))) // 1 + 1)
        max_x_abs = round(max(np.abs(list(qubit_x_pos))) // 1)

        self.flag_collection = []   # two rows of flag qubits
        for i in range(-max_x_abs,max_x_abs+1):
            for j in [-max_y_abs,max_y_abs]:
                qubit_new = Qubit(id_temp,(i,j))
                self.qubit_collection.append(qubit_new)
                self.qubit_network[qubit_new.pos] = qubit_new
                self.flag_network[(i,j//max_y_abs)] = qubit_new




    def _SE_compile(self, code_name: str) -> SE_instruction_set:
        code_temp = self.code_name_dict[code_name]
        z_ctrl_dict: dict[int,list] = {}
        z_targ_dict: dict[int,list] = {}
        x_ctrl_dict: dict[int,list] = {}
        x_targ_dict: dict[int,list] = {}
        mid_cycle_steps = []
        
        for qubit in code_temp.z_check_collection + code_temp.x_check_collection:
            for step_str in qubit.neighbor:
                if int(step_str) not in mid_cycle_steps:
                    mid_cycle_steps.append(int(step_str))
        mid_cycle_steps = sorted(mid_cycle_steps)

        for step in mid_cycle_steps:
            z_ctrl_dict[step] = []
            z_targ_dict[step] = []
            x_ctrl_dict[step] = []
            x_targ_dict[step] = []
            for qubit in code_temp.z_check_collection:
                if str(step) in qubit.neighbor:
                    ctrl_pos = qubit.neighbor[str(step)].pos
                    targ_pos = qubit.pos
                    z_ctrl_dict[step].append(self.qubit_network[ctrl_pos])
                    z_targ_dict[step].append(self.qubit_network[targ_pos])
            for qubit in code_temp.x_check_collection:
                if str(step) in qubit.neighbor:
                    ctrl_pos = qubit.pos
                    targ_pos = qubit.neighbor[str(step)].pos
                    x_ctrl_dict[step].append(self.qubit_network[ctrl_pos])
                    x_targ_dict[step].append(self.qubit_network[targ_pos])
        
        ret_se_instruction = SE_instruction_set(code_name, mid_cycle_steps,
                                                z_ctrl_dict, z_targ_dict,
                                                x_ctrl_dict, x_targ_dict)
        return ret_se_instruction


    def SE_round(self, code_name: str, cleaness: bool = False, post: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(z_checks_in_circ,noise=MR_noise)
        self.circ.reset_X(x_checks_in_circ,noise=MX_noise)
        self._add_to_active(z_checks_in_circ)
        self._add_to_active(x_checks_in_circ)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                                 noise=single_qubit_noise)
        # self.circ.single_qubit_gate('I', self.qubit_collection,
        #                             noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list: list[Qubit] = []
            targ_list: list[Qubit] = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            for qubit in ctrl_list + targ_list:
                counter = 0
                for another in ctrl_list + targ_list:
                    if qubit == another:
                        counter += 1
                if counter > 1:
                    print(qubit.pos,qubit.circuit_id)
            self.used_qubits = ctrl_list + targ_list
            self._idle_depo()
            
            self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ, 
                                noise = MR_noise)
        self.circ.measure_X(x_checks_in_circ, noise=MX_noise)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(x_checks_in_circ)
        self._pop_from_active(z_checks_in_circ)
        if post == True:
            for check in x_checks_in_circ + z_checks_in_circ:
                self.circ.detector([check], self.circ.time_temp - 1,
                                [check], self.circ.time_temp,
                            (check.pos[0],check.pos[1], self.circ.time_temp,self.stage))
        else:
            for check in x_checks_in_circ + z_checks_in_circ:
                self.circ.detector([check], self.circ.time_temp - 1,
                                [check], self.circ.time_temp,
                            (check.pos[0],check.pos[1], self.circ.time_temp))

        if post == True:
            self.stage_list.append(self.stage)
            self.stage += 1

        self.circ.clock_plus_one()
        self.circ.tick()
        
        pass

    def SE_round_z_corrected(self, code_name: str, cleaness: bool = False, post: bool = False,
                 z_sign_corrected: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR', self.noise_probability]
            Z_error = ['Z_ERROR', self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(z_checks_in_circ,noise=X_error)
        self.circ.reset_X(x_checks_in_circ,noise=Z_error)
        self._add_to_active(z_checks_in_circ)
        self._add_to_active(x_checks_in_circ)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list: list[Qubit] = []
            targ_list: list[Qubit] = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=tq_noise)
            for qubit in ctrl_list + targ_list:
                counter = 0
                for another in ctrl_list + targ_list:
                    if qubit == another:
                        counter += 1
                if counter > 1:
                    print(qubit.pos,qubit.circuit_id)
            self.used_qubits = ctrl_list + targ_list
            self._idle_depo(cleaness)
            
            self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ, 
                                noise = X_error)
        self.circ.measure_X(x_checks_in_circ, noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        self._pop_from_active(x_checks_in_circ)
        self._pop_from_active(z_checks_in_circ)
        if z_sign_corrected == False:
            if post == True:
                for check in x_checks_in_circ + z_checks_in_circ:
                    self.circ.detector([check], self.circ.time_temp - 1,
                                    [check], self.circ.time_temp,
                                (check.pos[0],check.pos[1], self.circ.time_temp,self.stage))
            else:
                for check in x_checks_in_circ + z_checks_in_circ:
                    self.circ.detector([check], self.circ.time_temp - 1,
                                    [check], self.circ.time_temp,
                                (check.pos[0],check.pos[1], self.circ.time_temp))
        elif z_sign_corrected == True:
            if post == True:
                for check in x_checks_in_circ:
                    self.circ.detector([check], self.circ.time_temp - 1,
                                    [check], self.circ.time_temp,
                                (check.pos[0],check.pos[1], self.circ.time_temp,self.stage))
                for check in z_checks_in_circ:
                    self.circ.detector_qt_list([check],[self.circ.time_temp],
                                (check.pos[0],check.pos[1],self.circ.time_temp,self.stage))
            else:
                for check in x_checks_in_circ:
                    self.circ.detector([check], self.circ.time_temp - 1,
                                    [check], self.circ.time_temp,
                                (check.pos[0],check.pos[1], self.circ.time_temp))
                for check in z_checks_in_circ:
                    self.circ.detector_qt_list([check],[self.circ.time_temp],
                                (check.pos[0],check.pos[1],self.circ.time_temp))


        if post == True:
            self.stage_list.append(self.stage)
            self.stage += 1

        self.circ.clock_plus_one()
        self.circ.tick()
        
        pass


    def SE_round_z_corrected_rp3(self, code_name: str, cleaness: bool = False, post: bool = False,
                 z_sign_corrected: bool = True,
                 detector_val: bool = True) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR', self.noise_probability]
            Z_error = ['Z_ERROR', self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(z_checks_in_circ,noise=X_error)
        self.circ.reset_X(x_checks_in_circ,noise=Z_error)
        self._add_to_active(z_checks_in_circ)
        self._add_to_active(x_checks_in_circ)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                                 noise=single_qubit_noise)
        # self.circ.single_qubit_gate('I', self.qubit_collection,
        #                             noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list: list[Qubit] = []
            targ_list: list[Qubit] = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=tq_noise)
            for qubit in ctrl_list + targ_list:
                counter = 0
                for another in ctrl_list + targ_list:
                    if qubit == another:
                        counter += 1
                if counter > 1:
                    print(qubit.pos,qubit.circuit_id)
            self.used_qubits = ctrl_list + targ_list
            self._idle_depo(cleaness)
            
            self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ, 
                                noise = X_error)
        self.circ.measure_X(x_checks_in_circ, noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        self._pop_from_active(x_checks_in_circ)
        self._pop_from_active(z_checks_in_circ)
        if detector_val == False:
            if z_sign_corrected == False:
                if post == True:
                    for check in x_checks_in_circ + z_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,self.stage))
                else:
                    for check in x_checks_in_circ + z_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp))
            elif z_sign_corrected == True:
                if post == True:
                    for check in x_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,self.stage))
                    for check in z_checks_in_circ:
                        self.circ.detector_qt_list([check],[self.circ.time_temp],
                                    (check.pos[0],check.pos[1],self.circ.time_temp,self.stage))
                else:
                    for check in x_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp))
                    for check in z_checks_in_circ:
                        self.circ.detector_qt_list([check],[self.circ.time_temp],
                                    (check.pos[0],check.pos[1],self.circ.time_temp))
        else:
            if z_sign_corrected == False:
                if post == True:
                    for check in x_checks_in_circ + z_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,self.stage,0))
                else:
                    for check in x_checks_in_circ + z_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,0))
            elif z_sign_corrected == True:
                if post == True:
                    for check in x_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,self.stage,0))
                    for check in z_checks_in_circ:
                        if check.pos == (1.5,-1.5):
                            self.circ.detector_qt_list([check],[self.circ.time_temp],
                                        (check.pos[0],check.pos[1],self.circ.time_temp,self.stage,0))
                        else:
                            self.circ.detector_qt_list([check],[self.circ.time_temp],
                                        (check.pos[0],check.pos[1],self.circ.time_temp,self.stage,1))
                else:
                    for check in x_checks_in_circ:
                        self.circ.detector([check], self.circ.time_temp - 1,
                                        [check], self.circ.time_temp,
                                    (check.pos[0],check.pos[1], self.circ.time_temp,0))
                    for check in z_checks_in_circ:
                        if check.pos == (1.5,-1.5):
                            self.circ.detector_qt_list([check],[self.circ.time_temp],
                                        (check.pos[0],check.pos[1],self.circ.time_temp,0))
                        else:
                            self.circ.detector_qt_list([check],[self.circ.time_temp],
                                        (check.pos[0],check.pos[1],self.circ.time_temp,1))


        if post == True:
            self.stage_list.append(self.stage)
            self.stage += 1

        self.circ.clock_plus_one()
        self.circ.tick()
        
        pass


    def SE_round_circ(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(z_checks_in_circ,noise=MR_noise)
        self.circ.reset_X(x_checks_in_circ,noise=MX_noise)
        self._add_to_active(z_checks_in_circ+x_checks_in_circ)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        # self.circ.single_qubit_gate('H', x_checks_in_circ)
        # self.circ.single_qubit_gate('I', self.qubit_collection,
        #                             noise = single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            self.used_qubits = ctrl_list+targ_list
            self._idle_depo()
            self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ,noise=MR_noise)
        self.circ.measure_X(x_checks_in_circ,noise=MX_noise)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()


    def X_init(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
            R_noise = ['X_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
            R_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        data_qubits_in_circ: list[Qubit] = []
        z_checks_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        
        self.circ.reset(data_qubits_in_circ + z_checks_in_circ \
                        + x_checks_in_circ, 
                        noise = R_noise)
        self.circ.single_qubit_gate('H', x_checks_in_circ + data_qubits_in_circ,
                                         noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)
        for qubit in x_checks_in_circ:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass

    def Z_init(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
            R_noise = ['X_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
            R_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        data_qubits_in_circ: list[Qubit] = []
        z_checks_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        
        self.circ.reset(data_qubits_in_circ + z_checks_in_circ \
                        + x_checks_in_circ, 
                        noise = MR_noise)
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                         noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)
        for qubit in z_checks_in_circ:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass

    def X_meas(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            MZ_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            MZ_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        data_qubits_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])

        self.circ.single_qubit_gate('H',data_qubits_in_circ,
                                    noise= single_qubit_noise)
        self.circ.measure(data_qubits_in_circ,
                          noise = MZ_noise)
        for x_check in code_temp.x_check_collection:
            check_pos = x_check.pos
            x_check_in_circ = self.qubit_network[check_pos]
            neighbors_in_circ = []
            for neighbor in x_check.neighbor.values():
                neighbors_in_circ.append(self.qubit_network[neighbor.pos])
            self.circ.detector([x_check_in_circ], self.circ.time_temp - 1,
                        neighbors_in_circ,
                        self.circ.time_temp,
                        (x_check_in_circ.pos[0],x_check_in_circ.pos[1], self.circ.time_temp))
        
        logic_x_in_circ = []
        for qubit in code_temp.logic_x_collection:
            logic_x_in_circ.append(self.qubit_network[qubit.pos])
        self.circ.observable(logic_x_in_circ,self.circ.time_temp)

        pass

    def Z_meas(self, code_name: str) -> None:
        code_temp = self.code_name_dict[code_name]
        data_qubits_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])

        self.circ.measure(data_qubits_in_circ,
                          noise = ['X_ERROR', self.noise_probability])
        for z_check in code_temp.z_check_collection:
            check_pos = z_check.pos
            z_check_in_circ = self.qubit_network[check_pos]
            neighbors_in_circ = []
            for neighbor in z_check.neighbor.values():
                neighbors_in_circ.append(self.qubit_network[neighbor.pos])
            self.circ.detector([z_check_in_circ], self.circ.time_temp - 1,
                        neighbors_in_circ,
                        self.circ.time_temp,
                        (z_check_in_circ.pos[0],z_check_in_circ.pos[1], self.circ.time_temp))
        
        logic_z_in_circ = []
        for qubit in code_temp.logic_z_collection:
            logic_z_in_circ.append(self.qubit_network[qubit.pos])
        self.circ.observable(logic_z_in_circ,self.circ.time_temp)
        pass


    def check_circuit_distance(self) -> int:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        explained_error = circuit.search_for_undetectable_logical_errors(
                            dont_explore_detection_event_sets_with_size_above=6,
                            dont_explore_edges_with_degree_above=3,
                            dont_explore_edges_increasing_symptom_degree=False)
        print(len(explained_error))
        for item in explained_error:
            print(item)
        return len(explained_error)
    
    def check_circuit_distance_graphlike(self) -> int:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        explained_error = circuit.shortest_graphlike_error()
        print(len(explained_error))
        for item in explained_error:
            print(item)
        return len(explained_error)
    

        
    def get_stim_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        return circuit
    

    def rp_to_sc(self, rp_name: str, sc_name: str) -> None:
        '''
        One-step growth from a smaller RP-surface code patch to 
        a larger rotated surface code patch.
        Growth: middle-out 
        '''
        tq_noise = ['DEPOLARIZE2', self.noise_probability]
        X_error = ['X_ERROR',self.noise_probability]
        Z_error = ['Z_ERROR',self.noise_probability]
        rp_code = self.code_name_dict[rp_name]
        sc_code = self.code_name_dict[sc_name]
        d_rp = rp_code.distance
        d_sc = sc_code.distance
        if d_rp >= d_sc:
            raise ValueError('rp code should be smaller')
        # step 1 reset 
        R_qubit_list = []
        RX_qubit_list = []
        for y in range(d_rp//2+1,d_sc//2+1):
            for x in range(-y,y):
                R_qubit_list.append(self.qubit_network[(x,-y)])
                R_qubit_list.append(self.qubit_network[(-x,y)])
        for x in range(d_rp//2+1,d_sc//2+1):
            for y in range(-x+1,x+1):
                RX_qubit_list.append(self.qubit_network[(-x,y)])
                RX_qubit_list.append(self.qubit_network[(x,-y)])

        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(R_qubit_list + z_checks_in_circ,
                        noise = ['X_ERROR', self.noise_probability])
        self.circ.reset_X(RX_qubit_list + x_checks_in_circ,
                          noise = ['Z_ERROR', self.noise_probability])
        self._add_to_active(R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ)
        self._idle_depo()
        self.circ.tick()


        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list

        for i in se_instruction.mid_cycle_steps:

            self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                        targ_dict[i],
                                        noise=tq_noise)
            self.used_qubits = ctrl_dict[i]+targ_dict[i]
            self._idle_depo()
            self.circ.tick()

        self.circ.measure(z_checks_in_circ,noise=X_error)
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()
        
        # detector specification:
        for check in sc_code.x_check_collection + sc_code.z_check_collection:
            qubit = self.qubit_network[check.pos]
            if np.abs(qubit.pos[0]) < d_rp // 2 and \
               np.abs(qubit.pos[1]) < d_rp // 2:
                self.circ.detector_qt_list([qubit, qubit],
                                [self.circ.time_temp-1, 
                                 self.circ.time_temp],
                         (qubit.pos[0],qubit.pos[1],self.circ.time_temp))
            elif np.abs(qubit.pos[0]) < d_rp // 2 and \
                 qubit.pos[1] == -(d_rp//2) - 0.5 and \
                 (qubit.pos[0]-0.5 + d_rp//2) % 2 == 0:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit,qubit_inv],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) < d_rp // 2 and \
                 qubit.pos[0] == -(d_rp//2) - 0.5 and \
                 (qubit.pos[1]+0.5-d_rp//2) % 2 == 0:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit,qubit_inv],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) > d_rp//2 + 1 and \
                 np.abs(qubit.pos[1]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[1]) > np.abs(qubit.pos[0]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) > d_rp//2 + 1 and \
                 np.abs(qubit.pos[0]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[0]) > np.abs(qubit.pos[1]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))            
        self.circ.clock_plus_one()
        self.circ.tick()        
        
        pass

    def color_to_rp3_to_sc_bell(self, rp_name: str, sc_name: str) -> None:
        

        tq_noise = ['DEPOLARIZE2', self.noise_probability]
        X_error = ['X_ERROR',self.noise_probability]
        Z_error = ['Z_ERROR',self.noise_probability]


        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
      
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1]
        self._idle_depo()
        self.circ.tick()

        rp_code = self.code_name_dict[rp_name]
        sc_code = self.code_name_dict[sc_name]
        d_rp = rp_code.distance
        d_sc = sc_code.distance
        if d_rp >= d_sc:
            raise ValueError('rp code should be smaller')
        if d_sc - d_rp < 4:
            raise ValueError('not supported. growth too small')
        # step 1 reset 
        R_qubit_list = []
        RX_qubit_list = []
        bell_c_list = []
        bell_t_list = []
        y_1_larger = d_rp//2 +1
        for x in range(-y_1_larger, y_1_larger):
            bell_c_list.append(self.qubit_network[(x,-y_1_larger)])
            bell_t_list.append(self.qubit_network[(-x,y_1_larger)])
        for y in range(-y_1_larger+1,y_1_larger+1):
            bell_c_list.append(self.qubit_network[(-y_1_larger,y)])
            bell_t_list.append(self.qubit_network[(y_1_larger,-y)])
        for y in range(y_1_larger+1,d_sc//2+1):
            for x in range(-y,y):
                R_qubit_list.append(self.qubit_network[(x,-y)])
                R_qubit_list.append(self.qubit_network[(-x,y)])
        for x in range(y_1_larger+1,d_sc//2+1):
            for y in range(-x+1,x+1):
                RX_qubit_list.append(self.qubit_network[(-x,y)])
                RX_qubit_list.append(self.qubit_network[(x,-y)])

        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])


        self.circ.reset(R_qubit_list + z_checks_in_circ + bell_t_list,
                        noise = ['X_ERROR', self.noise_probability])
        self.circ.reset_X(RX_qubit_list + x_checks_in_circ + bell_c_list,
                          noise = ['Z_ERROR', self.noise_probability])
        self._add_to_active(R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ
                            +bell_t_list+bell_c_list)
        self.used_qubits = R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ \
                            +bell_t_list+bell_c_list
        self._idle_depo()
        self.circ.tick()


        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1]+bell_c_list,
                                 four_packs[0]+four_packs[3]+bell_t_list,
                                 noise=tq_noise)
        self.used_qubits = bell_c_list + bell_t_list + \
                            four_packs[0] + four_packs[1] + four_packs[2] + four_packs[3]
        self._idle_depo()
        self.circ.tick()
        

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list


        first_step = se_instruction.mid_cycle_steps[0]




        # layer 1: bell pair making
        self.circ.two_qubit_gate('CNOT',four_packs[2] + ctrl_dict[first_step],
                                 four_packs[3] + targ_dict[first_step],
                                 noise = tq_noise)
        self.used_qubits = ctrl_dict[first_step]+targ_dict[first_step] \
                            + four_packs[2] + four_packs[3]
        self._idle_depo()
        self.circ.tick()


        self.circ.measure_X(four_packs[2],Z_error)
        self.circ.measure(four_packs[3],X_error)
        self._pop_from_active(four_packs[2]+four_packs[3])
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo()
        self.circ.tick()


        for qubit in four_packs[2] + four_packs[3]:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                    (qubit.pos[0], qubit.pos[1],
                        self.circ.time_temp, self.stage))
        
        self.stage_list.append(self.stage)
        self.stage += 1

        

        

    

        for i in se_instruction.mid_cycle_steps:
            if i != first_step:
                self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                            targ_dict[i],
                                            noise=tq_noise)
                self.used_qubits = ctrl_dict[i]+targ_dict[i]
                self._idle_depo()
                self.circ.tick()

        self.circ.measure(z_checks_in_circ,noise=X_error)
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()

        # anc positions in rp code:
        anc_pos_set = set()
        for check in rp_code.x_check_collection + rp_code.z_check_collection:
            anc_pos_set.add(check.pos)
        
        # detector specification:
        for check in sc_code.x_check_collection + sc_code.z_check_collection:
            qubit = self.qubit_network[check.pos]
            if np.abs(qubit.pos[0]) < d_rp // 2 and \
               np.abs(qubit.pos[1]) < d_rp // 2:
                self.circ.detector_qt_list([qubit, qubit],
                                [self.circ.time_temp-1, 
                                 self.circ.time_temp],
                         (qubit.pos[0],qubit.pos[1],self.circ.time_temp))
            elif qubit.pos in anc_pos_set:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit,qubit_inv],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) < (d_rp // 2) + 1 and \
                 qubit.pos[1] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) < (d_rp // 2) + 1 and \
                 qubit.pos[0] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[1]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[1]) > np.abs(qubit.pos[0]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[0]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[0]) > np.abs(qubit.pos[1]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))            
        self.circ.clock_plus_one()
        self.circ.tick()        
        

        pass

    def rp_to_sc_bell(self, rp_name: str, sc_name: str) -> None:
        '''
        One-step growth from a smaller RP-surface code patch to 
        a larger rotated surface code patch.
        Growth: middle-out 
        '''
        tq_noise = ['DEPOLARIZE2', self.noise_probability]
        X_error = ['X_ERROR',self.noise_probability]
        Z_error = ['Z_ERROR',self.noise_probability]
        rp_code = self.code_name_dict[rp_name]
        sc_code = self.code_name_dict[sc_name]
        d_rp = rp_code.distance
        d_sc = sc_code.distance
        if d_rp >= d_sc:
            raise ValueError('rp code should be smaller')
        if d_sc - d_rp < 4:
            raise ValueError('not supported. growth too small')
        # step 1 reset 
        R_qubit_list = []
        RX_qubit_list = []
        bell_c_list = []
        bell_t_list = []
        y_1_larger = d_rp//2 +1
        for x in range(-y_1_larger, y_1_larger):
            bell_c_list.append(self.qubit_network[(x,-y_1_larger)])
            bell_t_list.append(self.qubit_network[(-x,y_1_larger)])
        for y in range(-y_1_larger+1,y_1_larger+1):
            bell_c_list.append(self.qubit_network[(-y_1_larger,y)])
            bell_t_list.append(self.qubit_network[(y_1_larger,-y)])
        for y in range(y_1_larger+1,d_sc//2+1):
            for x in range(-y,y):
                R_qubit_list.append(self.qubit_network[(x,-y)])
                R_qubit_list.append(self.qubit_network[(-x,y)])
        for x in range(y_1_larger+1,d_sc//2+1):
            for y in range(-x+1,x+1):
                RX_qubit_list.append(self.qubit_network[(-x,y)])
                RX_qubit_list.append(self.qubit_network[(x,-y)])

        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(R_qubit_list + z_checks_in_circ + bell_t_list,
                        noise = ['X_ERROR', self.noise_probability])
        self.circ.reset_X(RX_qubit_list + x_checks_in_circ + bell_c_list,
                          noise = ['Z_ERROR', self.noise_probability])
        self._add_to_active(R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ
                            +bell_t_list+bell_c_list)
        self.used_qubits = R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ \
                            +bell_t_list+bell_c_list
        self._idle_depo()
        self.circ.tick()

        self.circ.two_qubit_gate('CNOT',bell_c_list,bell_t_list,noise=tq_noise)
        self.used_qubits = bell_c_list + bell_t_list
        self._idle_depo()
        self.circ.tick()

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list

        for i in se_instruction.mid_cycle_steps:

            self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                        targ_dict[i],
                                        noise=tq_noise)
            self.used_qubits = ctrl_dict[i]+targ_dict[i]
            self._idle_depo()
            self.circ.tick()

        self.circ.measure(z_checks_in_circ,noise=X_error)
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()

        # anc positions in rp code:
        anc_pos_set = set()
        for check in rp_code.x_check_collection + rp_code.z_check_collection:
            anc_pos_set.add(check.pos)
        
        # detector specification:
        for check in sc_code.x_check_collection + sc_code.z_check_collection:
            qubit = self.qubit_network[check.pos]
            if np.abs(qubit.pos[0]) < d_rp // 2 and \
               np.abs(qubit.pos[1]) < d_rp // 2:
                self.circ.detector_qt_list([qubit, qubit],
                                [self.circ.time_temp-1, 
                                 self.circ.time_temp],
                         (qubit.pos[0],qubit.pos[1],self.circ.time_temp))
            elif qubit.pos in anc_pos_set:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit,qubit_inv],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) < (d_rp // 2) + 1 and \
                 qubit.pos[1] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) < (d_rp // 2) + 1 and \
                 qubit.pos[0] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[1]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[1]) > np.abs(qubit.pos[0]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[0]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[0]) > np.abs(qubit.pos[1]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))            
        self.circ.clock_plus_one()
        self.circ.tick()        
        
        pass



    def Y_to_rp3(self, rp3_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            X_error = ['X_ERROR', self.noise_probability]
            Z_error = ['Z_ERROR', self.noise_probability]
            sq_noise = ['DEPOLARIZE1', self.noise_probability]
            tq_noise = ['DEPOLARIZE2', self.noise_probability]
        else:
            X_error = None
            Z_error = None
            sq_noise = None
            tq_noise = None
        rp3_code = self.code_name_dict[rp3_name]
        bell_pair_c = []
        bell_pair_t = []
        for i in [-1,0,1]:
            bell_pair_c.append(self.qubit_network[(i,-1)])
            bell_pair_t.append(self.qubit_network[(-i,1)])
        bell_pair_c.append(self.qubit_network[(-1,0)])
        bell_pair_t.append(self.qubit_network[(1,0)])
        # step 1 initialize the central qubit and the bell qubtis
        self.circ.reset_X(bell_pair_c, Z_error)
        self.circ.reset_X([self.qubit_network[(0,0)]],Z_error)
        self.circ.reset(bell_pair_t, X_error)
        self.circ.tick()
        self._add_to_active(bell_pair_c+bell_pair_t+[self.qubit_network[(0,0)]])
        self.circ.two_qubit_gate('CNOT', bell_pair_c, bell_pair_t,
                                 noise = tq_noise)
        self.circ.single_qubit_gate('S', [self.qubit_network[(0,0)]],
                                    noise = sq_noise)
        self.used_qubit_ids = bell_pair_c + bell_pair_t + [self.qubit_network[(0,0)]]
        self._idle_depo()
        self.circ.tick()
        # step 2 one round of stab measurement
        self.SE_round_circ(rp3_name, cleaness)
        for check in rp3_code.x_check_collection + rp3_code.z_check_collection:
            qubit = self.qubit_network[check.pos]
            if np.abs(qubit.pos[0]) > 1 or \
               np.abs(qubit.pos[1]) > 1:
                self.circ.detector_qt_list([qubit], [self.circ.time_temp],
                            (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage))
            elif qubit.pos[1] == -0.5:
                qubit_inv = self.qubit_network[(-qubit.pos[0],
                                                -qubit.pos[1])]
                self.circ.detector_qt_list([qubit, qubit_inv],
                                    [self.circ.time_temp,self.circ.time_temp],
                                    (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage))
        self.circ.clock_plus_one()
        self.stage_list.append(self.stage)
        self.stage += 1
        self.circ.tick()
        pass

    def Y_to_rp3_sign_correction(self, rp3_name: str, cleaness: bool = False,
                                 detector_val: bool = False) -> None:
        if cleaness == False:
            X_error = ['X_ERROR', self.noise_probability]
            Z_error = ['Z_ERROR', self.noise_probability]
            sq_noise = ['DEPOLARIZE1', self.noise_probability]
            tq_noise = ['DEPOLARIZE2', self.noise_probability]
        else:
            X_error = None
            Z_error = None
            tq_noise = None
            sq_noise = None
        rp3_code = self.code_name_dict[rp3_name]
        bell_pair_c = []
        bell_pair_t = []
        for i in [-1,0,1]:
            bell_pair_c.append(self.qubit_network[(i,-1)])
            bell_pair_t.append(self.qubit_network[(-i,1)])
        bell_pair_c.append(self.qubit_network[(-1,0)])
        bell_pair_t.append(self.qubit_network[(1,0)])
        # step 1 initialize the central qubit and the bell qubtis
        self.circ.reset_X(bell_pair_c, Z_error)
        self.circ.reset_X([self.qubit_network[(0,0)]],Z_error)
        self.circ.reset(bell_pair_t, X_error)
        self.circ.tick()
        self._add_to_active(bell_pair_c+bell_pair_t+[self.qubit_network[(0,0)]])
        self.circ.two_qubit_gate('CNOT', bell_pair_c, bell_pair_t,
                                 noise = tq_noise)
        self.circ.single_qubit_gate('S', [self.qubit_network[(0,0)]],
                                    noise = sq_noise)
        self.used_qubit_ids = bell_pair_c + bell_pair_t + [self.qubit_network[(0,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()
        # step 2 one round of stab measurement

        se_instruction = self.code_name_SE_instruction[rp3_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in rp3_code.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in rp3_code.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(z_checks_in_circ,noise=X_error)
        self.circ.reset_X(x_checks_in_circ,noise=Z_error)
        self._add_to_active(z_checks_in_circ+x_checks_in_circ)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        # self.circ.single_qubit_gate('H', x_checks_in_circ)
        # self.circ.single_qubit_gate('I', self.qubit_collection,
        #                             noise = single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=tq_noise)
            self.used_qubits = ctrl_list+targ_list
            self._idle_depo(cleaness)
            self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ,noise=X_error,
                          not_qubits=[self.qubit_network[(-0.5,0.5)],
                                      self.qubit_network[(0.5,-0.5)]])
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)

        self.circ.single_qubit_gate('X',[self.qubit_network[(-1,-1)]])

        self.circ.feedforward_pauli('X', 
                    meas_qubit_list=[self.qubit_network[(-0.5,0.5)],
                                     self.qubit_network[(0.5,-0.5)]],
                    corr_qubit_list=[self.qubit_network[(-1,0)],
                                     self.qubit_network[(1,0)]],
                    time_or_stage_list=[(self.circ.time_temp,'time'),
                                        (self.circ.time_temp,'time')])
        
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()

        if detector_val == False:
            for check in rp3_code.x_check_collection + rp3_code.z_check_collection:
                qubit = self.qubit_network[check.pos]
                if np.abs(qubit.pos[0]) > 1 or \
                    np.abs(qubit.pos[1]) > 1:
                    self.circ.detector_qt_list([qubit], [self.circ.time_temp],
                                (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage))
                elif qubit.pos[1] == -0.5:
                    qubit_inv = self.qubit_network[(-qubit.pos[0],
                                                    -qubit.pos[1])]
                    self.circ.detector_qt_list([qubit, qubit_inv],
                                        [self.circ.time_temp,self.circ.time_temp],
                                        (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage))
        else:
            for check in rp3_code.x_check_collection + rp3_code.z_check_collection:
                qubit = self.qubit_network[check.pos]
                if np.abs(qubit.pos[0]) > 1 or \
                    np.abs(qubit.pos[1]) > 1:
                    self.circ.detector_qt_list([qubit], [self.circ.time_temp],
                                (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage,0))
                elif qubit.pos[1] == -0.5:
                    qubit_inv = self.qubit_network[(-qubit.pos[0],
                                                    -qubit.pos[1])]
                    self.circ.detector_qt_list([qubit, qubit_inv],
                                        [self.circ.time_temp,self.circ.time_temp],
                                        (qubit.pos[0], qubit.pos[1], self.circ.time_temp, self.stage,0))
        self.circ.clock_plus_one()
        self.stage_list.append(self.stage)
        self.stage += 1
        self.circ.tick()
        pass

    def magic_Y_measurement(self, code_name: str) -> None:
        code = self.code_name_dict[code_name]
        logical_qubit = Qubit(-1, (0,0))
        self.logic_qubits.append(logical_qubit)
        stab_supps = []
        stab_opers = []
        registers = []
        for x_check in code.x_check_collection:
            supp = []
            oper = []
            for neighbor in x_check.neighbor.values():
                if neighbor == None:
                    continue
                neighbor_qubit = self.qubit_network[neighbor.pos]
                supp.append(neighbor_qubit)
                oper.append('X')
            stab_supps.append(tuple(supp))
            stab_opers.append(tuple(oper))
            registers.append(self.qubit_network[x_check.pos])

        for z_check in code.z_check_collection:
            supp = []
            oper = []
            for neighbor in z_check.neighbor.values():
                if neighbor == None:
                    continue
                neighbor_qubit = self.qubit_network[neighbor.pos]
                supp.append(neighbor_qubit)
                oper.append('Z')
            stab_supps.append(tuple(supp))
            stab_opers.append(tuple(oper))
            registers.append(self.qubit_network[z_check.pos])
        
        logic_supp = []
        logic_oper = []
        for code_qubit, code_oper in zip(code.logic_y_supp,
                                         code.logic_y_oper):
            qubit = self.qubit_network[code_qubit.pos]
            logic_supp.append(qubit)
            logic_oper.append(code_oper)
        stab_supps.append(tuple(logic_supp))
        stab_opers.append(tuple(logic_oper))
        registers.append(logical_qubit)

        self.circ.MPP(stab_supps, stab_opers, registers)
        for qubit in registers:
            if qubit == logical_qubit:
                self.circ.observable([qubit],self.circ.time_temp)
            else:
                self.circ.detector_qt_list([qubit, qubit],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp],
                    (qubit.pos[0], qubit.pos[1], self.circ.time_temp))

        self.circ.clock_plus_one()
        pass

    def magic_Y_measurement_post_selection(self, code_name: str) -> None:
        code = self.code_name_dict[code_name]
        logical_qubit = Qubit(-1, (0,0))
        self.logic_qubits.append(logical_qubit)
        stab_supps = []
        stab_opers = []
        registers = []
        for x_check in code.x_check_collection:
            supp = []
            oper = []
            for neighbor in x_check.neighbor.values():
                if neighbor == None:
                    continue
                neighbor_qubit = self.qubit_network[neighbor.pos]
                supp.append(neighbor_qubit)
                oper.append('X')
            stab_supps.append(tuple(supp))
            stab_opers.append(tuple(oper))
            registers.append(self.qubit_network[x_check.pos])

        for z_check in code.z_check_collection:
            supp = []
            oper = []
            for neighbor in z_check.neighbor.values():
                if neighbor == None:
                    continue
                neighbor_qubit = self.qubit_network[neighbor.pos]
                supp.append(neighbor_qubit)
                oper.append('Z')
            stab_supps.append(tuple(supp))
            stab_opers.append(tuple(oper))
            registers.append(self.qubit_network[z_check.pos])
        
        logic_supp = []
        logic_oper = []
        for code_qubit, code_oper in zip(code.logic_y_supp,
                                         code.logic_y_oper):
            qubit = self.qubit_network[code_qubit.pos]
            logic_supp.append(qubit)
            logic_oper.append(code_oper)
        stab_supps.append(tuple(logic_supp))
        stab_opers.append(tuple(logic_oper))
        registers.append(logical_qubit)

        self.circ.MPP(stab_supps, stab_opers, registers)
        for qubit in registers:
            if qubit == logical_qubit:
                self.circ.observable([qubit],self.circ.time_temp)
            else:
                self.circ.detector_qt_list([qubit, qubit],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp],
                    (qubit.pos[0], qubit.pos[1], self.circ.time_temp,
                     self.stage))
        
        self.circ.clock_plus_one()
        pass



    def rp3_to_rp5(self, rp3_name: str, rp5_name) -> None:
        '''
        One-step growth from a smaller RP-surface code patch to a 
        larger RP-surface code patch. 
        Restricted to the following case:
        d = 3 for the smaller patch and d = 5 for the larger patch.
        Bell pairs are used.
        '''
        tq_noise = ['DEPOLARIZE2',self.noise_probability]
        X_error = ['X_ERROR',self.noise_probability]
        Z_error = ['Z_ERROR',self.noise_probability]
        rp3_code = self.code_name_dict[rp3_name]
        rp5_code = self.code_name_dict[rp5_name]
        bell_pair_c = []
        bell_pair_t = []
        for i in [-2,-1,0,1,2]:
            bell_pair_c.append(self.qubit_network[(i,-2)])
            bell_pair_t.append(self.qubit_network[(-i,2)])
        for j in [-1,0,1]:
            bell_pair_c.append(self.qubit_network[(2,j)])
            bell_pair_t.append(self.qubit_network[(-2,-j)])
        
        code_temp = self.code_name_dict[rp5_name]
        se_instruction = self.code_name_SE_instruction[rp5_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        # step 1 initiate bell pairs
        self.circ.reset_X(bell_pair_c + x_checks_in_circ, 
                          noise=['Z_ERROR',self.noise_probability])
        self.circ.reset(bell_pair_t + z_checks_in_circ,
                        noise= ['X_ERROR', self.noise_probability])
        self._add_to_active(bell_pair_c+bell_pair_t+x_checks_in_circ+z_checks_in_circ)
        self.used_qubits = bell_pair_c + bell_pair_t + x_checks_in_circ + z_checks_in_circ
        self._idle_depo()
        self.circ.tick()

        self.circ.two_qubit_gate('CNOT', bell_pair_c, bell_pair_t,
                    noise = ['DEPOLARIZE2', self.noise_probability])
        self.used_qubits = bell_pair_c + bell_pair_t
        self._idle_depo()
        self.circ.tick()

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list

        for i in se_instruction.mid_cycle_steps:

            self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                        targ_dict[i],
                                        noise=tq_noise)
            self.used_qubits = ctrl_dict[i]+targ_dict[i]
            self._idle_depo()
            self.circ.tick()

        self.circ.measure(z_checks_in_circ,noise=X_error)
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()

        for qubit in rp5_code.x_check_collection + rp5_code.z_check_collection:
            qubit_circ = self.qubit_network[qubit.pos]
            if np.abs(qubit_circ.pos[0]) < 1 and \
               np.abs(qubit_circ.pos[1]) < 1:
                self.circ.detector(
                            [qubit_circ], self.circ.time_temp - 1,
                            [qubit_circ], self.circ.time_temp,
                            (qubit_circ.pos[0],qubit_circ.pos[1],
                             self.circ.time_temp, self.stage))
            elif qubit_circ.pos[1] == -1.5 and \
                    np.abs(qubit_circ.pos[0]) < 2:
                qubit_inv_pos = (-qubit_circ.pos[0],1.5)
                qubit_inv_circ = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list(
                            [qubit_circ, qubit_circ, qubit_inv_circ],
                            [self.circ.time_temp-1, self.circ.time_temp,
                             self.circ.time_temp],
                             (qubit_circ.pos[0], qubit_circ.pos[1],
                              self.circ.time_temp, self.stage))
            elif qubit_circ.pos[0] == -1.5 and \
                    np.abs(qubit_circ.pos[1]) < 1:
                qubit_inv_pos = (1.5,-qubit_circ.pos[1])
                qubit_inv_circ = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list(
                            [qubit_circ, qubit_circ, qubit_inv_circ],
                            [self.circ.time_temp-1, self.circ.time_temp,
                             self.circ.time_temp],
                             (qubit_circ.pos[0], qubit_circ.pos[1],
                              self.circ.time_temp, self.stage))
            elif np.abs(qubit_circ.pos[0]) > 2 or \
                 np.abs(qubit_circ.pos[1]) > 2:
                self.circ.detector_qt_list(
                            [qubit_circ], [self.circ.time_temp],
                            (qubit_circ.pos[0], qubit_circ.pos[1],
                             self.circ.time_temp, self.stage))
        



        self.circ.clock_plus_one()
        self.stage_list.append(self.stage)
        self.stage += 1
        self.circ.tick()
        pass


    def rp3_to_color_compact(self, rp3_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
        self.circ.reset_X(four_packs[2],Z_error)
        self.circ.reset(four_packs[3],X_error)
        self._add_to_active(four_packs[2]+four_packs[3])
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 1: bell pair making
        self.circ.two_qubit_gate('CNOT',four_packs[2],four_packs[3],
                                 noise = tq_noise)
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1],
                                 four_packs[0]+four_packs[3],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1] + four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()
        
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1]
        self._idle_depo(cleaness)
        self.circ.tick()

    def color_to_rp3_compact(self, rp3_name: str, cleaness: bool = False,
                             detector_val: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
        # self._idle_depo()
        
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1],
                                 four_packs[0]+four_packs[3],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1] + four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()
        
        # layer 1: bell pair making
        self.circ.two_qubit_gate('CNOT',four_packs[2],four_packs[3],
                                 noise = tq_noise)
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        self.circ.flag_measure_X(four_packs[2],self.stage,Z_error)
        self.circ.flag_measure(four_packs[3],self.stage,X_error)
        self._pop_from_active(four_packs[2]+four_packs[3])
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        if detector_val == False:
            for qubit in four_packs[2] + four_packs[3]:
                self.circ.flag_detector(qubit,self.stage, 
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage))
        else:
            for qubit in four_packs[2] + four_packs[3]:
                self.circ.flag_detector(qubit,self.stage, 
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage,0))
        
        self.stage_list.append(self.stage)
        self.stage += 1


    def color_to_rp3_compact_no_growth(self, rp3_name: str, 
                                       sc_name: str,
                                       cleaness: bool = False,
                             detector_val: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
        # self._idle_depo()
        
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1]
        self._idle_depo(cleaness)
        self.circ.tick()

        # reset sc code qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        


        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1],
                                 four_packs[0]+four_packs[3],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1] + four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()


        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list


        first_step = se_instruction.mid_cycle_steps[0]
        qubits_first_step = ctrl_dict[first_step] + targ_dict[first_step]

        
        # layer 1: bell pair making
        self.circ.two_qubit_gate('CNOT',four_packs[2],four_packs[3],
                                 noise = tq_noise)
        self.used_qubits = four_packs[2] + four_packs[3] + qubits_first_step
        self._idle_depo(cleaness)
        self.circ.tick()

        self.circ.flag_measure_X(four_packs[2],self.stage,Z_error)
        self.circ.flag_measure(four_packs[3],self.stage,X_error)
        self._pop_from_active(four_packs[2]+four_packs[3])
        self.used_qubits = four_packs[2] + four_packs[3] + qubits_first_step
        self._idle_depo(cleaness)
        self.circ.tick()

        if detector_val == False:
            for qubit in four_packs[2] + four_packs[3]:
                self.circ.flag_detector(qubit,self.stage, 
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage))
        else:
            for qubit in four_packs[2] + four_packs[3]:
                self.circ.flag_detector(qubit,self.stage, 
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage,0))
        
        self.stage_list.append(self.stage)
        self.stage += 1
    
    
    def rp3_to_color_and_MPP(self, detector_val: bool = True) -> None:
        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
        self.circ.reset_X(four_packs[2])
        self.circ.reset(four_packs[3])
        self.circ.tick()

        # layer 1: bell pair making
        self.circ.two_qubit_gate('CNOT',four_packs[2],four_packs[3])
        self.circ.tick()

        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1],
                                 four_packs[0]+four_packs[3])
        self.circ.tick()
        
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],four_packs[1])
        self.circ.tick()

        meas_supps: list[tuple[Qubit]] = []
        meas_opers: list[tuple[str]] = []
        qubit_register: list[Qubit] = []
        # face stabilizers
        for i in [-1,0,1]:
            meas_supps_temp = [self.qubit_network[(i,-1)],
                               self.qubit_network[(i,1)],
                               self.flag_network[(i,-1)],
                               self.flag_network[(i,1)]]
            meas_opers_temp = ['X','X','X','X']

            meas_supps.append(tuple(meas_supps_temp))
            meas_opers.append(tuple(meas_opers_temp))
            qubit_register.append(self.flag_network[(i,-1)])

            meas_opers_temp = ['Z','Z','Z','Z']
            meas_supps.append(tuple(meas_supps_temp))
            meas_opers.append(tuple(meas_opers_temp))
            qubit_register.append(self.flag_network[(i,1)])
        # rp surface code stabs
        meas_supps_temp = [self.qubit_network[(-1,0)],
                           self.qubit_network[(0,0)],
                           self.qubit_network[(-1,-1)],
                           self.qubit_network[(-1,1)],
                           self.qubit_network[(0,-1)],
                           self.qubit_network[(0,1)]]
        meas_opers_temp = ['X' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-0.5,-0.5)])

        meas_opers_temp = ['Z' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-0.5,0.5)])

        meas_supps_temp = [self.qubit_network[(0,0)],
                           self.qubit_network[(1,0)],
                           self.qubit_network[(0,1)],
                           self.flag_network[(0,1)],
                           self.qubit_network[(1,1)],
                           self.flag_network[(1,1)]]
        meas_opers_temp = ['X' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(0.5,0.5)])

        meas_opers_temp = ['Z' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(0.5,-0.5)])


        meas_supps_temp = [self.qubit_network[(-1,1)],
                           self.flag_network[(-1,1)],
                           self.qubit_network[(1,-1)],
                           self.qubit_network[(1,1)],
                           self.qubit_network[(0,-1)],
                           self.flag_network[(0,1)]]
        meas_opers_temp = ['X' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(0.5,-1.5)])

        meas_opers_temp = ['Z' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-0.5,-1.5)])

        meas_supps_temp = [self.qubit_network[(-1,0)],
                           self.qubit_network[(1,0)],
                           self.qubit_network[(-1,1)],
                           self.flag_network[(-1,1)],
                           self.qubit_network[(1,-1)],
                           self.qubit_network[(1,1)]]
        meas_opers_temp = ['X' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-1.5,0.5)])

        meas_opers_temp = ['Z' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-1.5,-0.5)])


        meas_supps_temp = [self.qubit_network[(-1,-1)],
                           self.qubit_network[(-1,1)],
                           self.qubit_network[(1,1)],
                           self.flag_network[(1,1)]]
        meas_opers_temp = ['X' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(-1.5,-1.5)])

        meas_opers_temp = ['Z' for _ in meas_supps_temp]
        meas_supps.append(tuple(meas_supps_temp))
        meas_opers.append(tuple(meas_opers_temp))
        qubit_register.append(self.qubit_network[(1.5,-1.5)])

        meas_supps_temp = []
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                meas_supps_temp.append(self.qubit_network[(i,j)])
        for i in [-1,0,1]:
            for j in [-1,1]:
                meas_supps_temp.append(self.flag_network[(i,j)])
        
        meas_opers_temp = ['Y' for _ in meas_supps_temp]
        meas_supps.append(meas_supps_temp)
        meas_opers.append(meas_opers_temp)
        logical_qubit = Qubit(-1, (0,0))
        self.logic_qubits.append(logical_qubit)
        qubit_register.append(logical_qubit)



        self.circ.MPP(meas_supps=meas_supps,
                      meas_opers=meas_opers,
                      qubit_register=qubit_register)
        

        if detector_val == True:
            for qubit in qubit_register:
                if qubit in self.flag_network.values():
                    self.circ.detector_qt_list([qubit],
                                            [self.circ.time_temp],
                        (qubit.pos[0], qubit.pos[1], self.circ.time_temp,0))
                elif qubit not in self.logic_qubits:
                    self.circ.detector_qt_list([qubit,qubit],
                                        [self.circ.time_temp-1,
                                        self.circ.time_temp],
                            (qubit.pos[0],qubit.pos[1],self.circ.time_temp,0))
                else:
                    self.circ.observable([logical_qubit],
                                        self.circ.time_temp)
                    pass


    
    def color_to_rp3_to_rp5(self, rp3_name: str, rp5_name, 
                            post: bool = True,
                            cleaness: bool = False) -> None:
        '''
        One-step growth from a smaller RP-surface code patch to a 
        larger RP-surface code patch. 
        Restricted to the following case:
        d = 3 for the smaller patch and d = 5 for the larger patch.
        Bell pairs are used.
        '''
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None
        
        rp3_code = self.code_name_dict[rp3_name]
        rp5_code = self.code_name_dict[rp5_name]

        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        
        code_temp = self.code_name_dict[rp5_name]
        se_instruction = self.code_name_SE_instruction[rp5_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        bell_pair_c = []
        bell_pair_t = []
        for i in [-2,-1,0,1,2]:
            bell_pair_c.append(self.qubit_network[(i,-2)])
            bell_pair_t.append(self.qubit_network[(-i,2)])
        for j in [-1,0,1]:
            bell_pair_c.append(self.qubit_network[(2,j)])
            bell_pair_t.append(self.qubit_network[(-2,-j)])

        
                
        # layer 3: stab inter-weaving
        self.circ.two_qubit_gate('CNOT',four_packs[0],
                                 four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1] 
        self._idle_depo(cleaness)
        self.circ.tick()


        # reset data for growth and ancillas 
        self.circ.reset(z_checks_in_circ + bell_pair_t,noise=X_error)
        self.circ.reset_X(x_checks_in_circ + bell_pair_c,noise=Z_error)
        self._add_to_active(z_checks_in_circ+x_checks_in_circ
                            +bell_pair_c+bell_pair_t)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ \
                            + bell_pair_t + bell_pair_c
        self._idle_depo(cleaness)
        self.circ.tick()


        # layer 2: stab expansion 
        self.circ.two_qubit_gate('CNOT',four_packs[2]+four_packs[1]+bell_pair_c,
                                 four_packs[0]+four_packs[3]+bell_pair_t,
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + four_packs[1] + four_packs[2] + four_packs[3] \
                            + bell_pair_c + bell_pair_t
        self._idle_depo(cleaness)
        self.circ.tick()

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list
        
        first_step = min(se_instruction.mid_cycle_steps)

        # layer 1: bell pair making and first step of SE.
        self.circ.two_qubit_gate('CNOT',four_packs[2] + ctrl_dict[first_step],
                                 four_packs[3] + targ_dict[first_step],
                                 noise = tq_noise)
        self.used_qubits = four_packs[2] + four_packs[3] \
                            + ctrl_dict[first_step] + targ_dict[first_step]
        self._idle_depo(cleaness)
        self.circ.tick()
        
        self.circ.flag_measure_X(four_packs[2],self.stage,Z_error)
        self.circ.flag_measure(four_packs[3],self.stage,X_error)
        self._pop_from_active(four_packs[2]+four_packs[3])
        self.used_qubits = four_packs[2] + four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        for qubit in four_packs[2] + four_packs[3]:
            self.circ.flag_detector(qubit, self.stage, 
                    (qubit.pos[0], qubit.pos[1],
                        self.circ.time_temp, self.stage))

        self.stage_list.append(self.stage)
        self.stage += 1
        



        for i in se_instruction.mid_cycle_steps:
            if i == first_step:
                continue
            else:

                self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                            targ_dict[i],
                                            noise=tq_noise)
                self.used_qubits = ctrl_dict[i]+targ_dict[i]
                self._idle_depo(cleaness)
                self.circ.tick()
        
        # self.circ.single_qubit_gate('H', x_checks_in_circ,
        #                             noise=single_qubit_noise)
        self.circ.measure(z_checks_in_circ,noise=X_error,
                    not_qubits=[self.qubit_network[(-0.5,-1.5)],
                                self.qubit_network[(1.5,0.5)]])
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.circ.single_qubit_gate('X',[self.qubit_network[(-1,-1)],
                                    self.qubit_network[(-2,1)],
                                    self.qubit_network[(-1,2)]])
        self.circ.feedforward_pauli('X',
                meas_qubit_list=[self.qubit_network[(-1.5,-0.5)],
                                 self.qubit_network[(1.5,0.5)],
                                 self.qubit_network[(-0.5,-1.5)],
                                 self.qubit_network[(0.5,1.5)],
                                 self.qubit_network[(-1.5,1.5)],
                                 self.qubit_network[(1.5,-1.5)]],
                corr_qubit_list=[self.qubit_network[(-2,0)],
                                 self.qubit_network[(2,0)],
                                 self.qubit_network[(0,-2)],
                                 self.qubit_network[(0,2)],
                                 self.qubit_network[(-2,2)],
                                 self.qubit_network[(2,-2)]],
                time_or_stage_list=[(self.circ.time_temp,'time')]*6)
        
        
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo(cleaness)
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()


        for qubit in rp5_code.x_check_collection + rp5_code.z_check_collection:
            qubit_circ = self.qubit_network[qubit.pos]
            if np.abs(qubit_circ.pos[0]) < 1 and \
               np.abs(qubit_circ.pos[1]) < 1:
                self.circ.detector(
                            [qubit_circ], self.circ.time_temp - 1,
                            [qubit_circ], self.circ.time_temp,
                            (qubit_circ.pos[0],qubit_circ.pos[1],
                             self.circ.time_temp, self.stage))
            elif qubit_circ.pos[1] == -1.5 and \
                    np.abs(qubit_circ.pos[0]) < 2:
                qubit_inv_pos = (-qubit_circ.pos[0],1.5)
                qubit_inv_circ = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list(
                            [qubit_circ, qubit_circ, qubit_inv_circ],
                            [self.circ.time_temp-1, self.circ.time_temp,
                             self.circ.time_temp],
                             (qubit_circ.pos[0], qubit_circ.pos[1],
                              self.circ.time_temp, self.stage))
            elif qubit_circ.pos[0] == -1.5 and \
                    np.abs(qubit_circ.pos[1]) < 1:
                qubit_inv_pos = (1.5,-qubit_circ.pos[1])
                qubit_inv_circ = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list(
                            [qubit_circ, qubit_circ, qubit_inv_circ],
                            [self.circ.time_temp-1, self.circ.time_temp,
                             self.circ.time_temp],
                             (qubit_circ.pos[0], qubit_circ.pos[1],
                              self.circ.time_temp, self.stage))
            elif np.abs(qubit_circ.pos[0]) > 2 or \
                 np.abs(qubit_circ.pos[1]) > 2:
                self.circ.detector_qt_list(
                            [qubit_circ], [self.circ.time_temp],
                            (qubit_circ.pos[0], qubit_circ.pos[1],
                             self.circ.time_temp, self.stage))
                


        self.circ.clock_plus_one()
        self.stage_list.append(self.stage)
        self.stage += 1
        self.circ.tick()
        pass

    def d3_color_T_check_compact(self, rp3_name: str, cleaness: bool = False,
                                 detector_val: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        four_packs = [[self.qubit_network[(i,-1)] for i in [-1,0,1]], 
                      [self.qubit_network[(i,1)] for i in [-1,0,1]],
                      [self.flag_network[(i,-1)] for i in [-1,0,1]],
                      [self.flag_network[(i,1)] for i in [-1,0,1]]]
        central_line_qubits = [self.qubit_network[(i,0)] for i in [-1,0,1]]
        all_srp_code_qubits = []
        all_srp_code_qubits += central_line_qubits
        for pack in four_packs:
            all_srp_code_qubits += pack

        shadow_qubits = []
        real_qubits = []
        for i in [-1,0,1]:
            for j in [-1,1]:
                shadow_qubits.append(self.qubit_network[(i-0.5,j-0.5)])
                real_qubits.append(self.qubit_network[(i,j)])
        for i in [-1,0,1]:
            shadow_qubits.append(self.qubit_network[(i-0.5,-0.5)])
            real_qubits.append(self.qubit_network[(i,0)])
        

        # S-gate
        self.circ.single_qubit_gate('S',all_srp_code_qubits,
                                    noise=sq_noise)
        self.used_qubits = all_srp_code_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        # initialize shadow qubits:
        self.circ.reset_X(shadow_qubits,Z_error)
        self._add_to_active(shadow_qubits)
        self.used_qubits = shadow_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        # bell pair + X retraction-1
        self.circ.two_qubit_gate('CNOT',shadow_qubits,
                                 real_qubits,
                                 noise=tq_noise)
        self.used_qubits = shadow_qubits \
                           + real_qubits 
        self._idle_depo(cleaness)
        self.circ.tick()

        # X retraction-1
        self.circ.two_qubit_gate('CNOT',four_packs[1],four_packs[3],
                                 noise=tq_noise)
        self.used_qubits = four_packs[1]+four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X retraction-2 
        self.circ.two_qubit_gate('CNOT', four_packs[0]+central_line_qubits,
                                 four_packs[2]+four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + central_line_qubits \
                           + four_packs[2] + four_packs[1]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X retraction-3
        self.circ.two_qubit_gate('CNOT', central_line_qubits,
                                 four_packs[0], noise= tq_noise)    
        self.used_qubits = central_line_qubits + four_packs[0]
        self._idle_depo(cleaness)
        self.circ.tick() 

        # X retract_central-1
        self.circ.two_qubit_gate('CNOT',[self.qubit_network[(0,0)]],
                                 [self.qubit_network[(-1,0)]],
                                 noise=tq_noise)
        self.used_qubits = [self.qubit_network[(0,0)],
                            self.qubit_network[(-1,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X retract_central-2 
        self.circ.two_qubit_gate('CNOT',[self.qubit_network[(0,0)]],
                                 [self.qubit_network[(1,0)]],
                                 noise=tq_noise)
        self.used_qubits = [self.qubit_network[(0,0)],
                            self.qubit_network[(1,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        # meas-X
        self.circ.flag_measure_X([self.qubit_network[(0,0)]],
                                 self.stage,noise=Z_error)
        
        self.used_qubits = [self.qubit_network[(0,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()
        if detector_val == False:
            self.circ.flag_detector(self.qubit_network[(0,0)], 
                                    self.stage, 
                        (0, 0,
                            self.circ.time_temp, self.stage))
        else:
            self.circ.flag_detector(self.qubit_network[(0,0)], 
                                    self.stage, 
                        (0, 0,
                            self.circ.time_temp, self.stage,0))
        self.stage += 1

        # reset-X
        self.circ.reset_X([self.qubit_network[(0,0)]],noise=Z_error)
        self.used_qubits = [self.qubit_network[(0,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X grow_central-1 
        self.circ.two_qubit_gate('CNOT',[self.qubit_network[(0,0)]],
                                 [self.qubit_network[(1,0)]],
                                 noise=tq_noise)
        self.used_qubits = [self.qubit_network[(0,0)],
                            self.qubit_network[(1,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X grow_central-2
        self.circ.two_qubit_gate('CNOT',[self.qubit_network[(0,0)]],
                                 [self.qubit_network[(-1,0)]],
                                 noise=tq_noise)
        self.used_qubits = [self.qubit_network[(0,0)],
                            self.qubit_network[(-1,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X growth-1
        self.circ.two_qubit_gate('CNOT', central_line_qubits,
                                 four_packs[0], noise= tq_noise)    
        self.used_qubits = central_line_qubits + four_packs[0]
        self._idle_depo(cleaness)
        self.circ.tick() 

        # X growth-2 
        self.circ.two_qubit_gate('CNOT', four_packs[0]+central_line_qubits,
                                 four_packs[2]+four_packs[1],
                                 noise=tq_noise)
        self.used_qubits = four_packs[0] + central_line_qubits \
                           + four_packs[2] + four_packs[1]
        self._idle_depo(cleaness)
        self.circ.tick()

        # X growth-3
        self.circ.two_qubit_gate('CNOT',four_packs[1],four_packs[3],
                                 noise=tq_noise)
        self.used_qubits = four_packs[1]+four_packs[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # bell pair unmaking 
        self.circ.two_qubit_gate('CNOT',shadow_qubits,
                                 real_qubits,
                                 noise=tq_noise)
        self.used_qubits = shadow_qubits \
                           + real_qubits 
        self._idle_depo(cleaness)
        self.circ.tick()

        # measure shadow qubits:
        self.circ.flag_measure_X(shadow_qubits,
                                 self.stage,Z_error)
        self._pop_from_active(shadow_qubits)
        self.used_qubits = shadow_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        # S-gate
        self.circ.single_qubit_gate('S_DAG',all_srp_code_qubits,
                                    noise=sq_noise)
        self.used_qubits = all_srp_code_qubits
        self._idle_depo(cleaness)
        self.circ.tick()


        
        if detector_val == False:
            for qubit in shadow_qubits:
                self.circ.flag_detector(qubit, self.stage,
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage))
        else:
            for qubit in shadow_qubits:
                self.circ.flag_detector(qubit, self.stage,
                        (qubit.pos[0], qubit.pos[1],
                            self.circ.time_temp, self.stage,0))


        self.stage_list.append(self.stage)
        self.stage += 1
        
        pass

        

    def rp5_to_color_new(self, cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        six_packs = []

        for i in [-2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-1)])
        six_pack.append(self.qubit_network[(-1,1)])
        six_pack.append(self.qubit_network[(0,-1)])
        six_pack.append(self.qubit_network[(0,1)])
        six_pack.append(self.flag_network[(-1,-1)])
        six_pack.append(self.flag_network[(-1,1)])
        six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-2)])
        six_pack.append(self.qubit_network[(-1,2)])
        six_pack.append(self.qubit_network[(0,-2)])
        six_pack.append(self.qubit_network[(0,2)])
        six_pack.append(self.flag_network[(0,-1)])
        six_pack.append(self.flag_network[(0,1)])
        six_packs.append(six_pack)
        for i in [1,2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        
        six_pack_col_list = []
        for i in range(6):
            col_temp = []
            for six_pack in six_packs:
                col_temp.append(six_pack[i])
            six_pack_col_list.append(col_temp)

        # layer 1: Bell pair 
        self.circ.reset_X(six_pack_col_list[4],
                            noise=Z_error)
        self.circ.reset(six_pack_col_list[5])
        self.circ.single_qubit_gate('X',six_pack_col_list[5],
                                    noise=X_error)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5]
        self._add_to_active(six_pack_col_list[4] + six_pack_col_list[5])
        self._idle_depo(cleaness)
        self.circ.tick()

        # bell pair making

        self.circ.two_qubit_gate('CNOT', 
                                 six_pack_col_list[4],
                                 six_pack_col_list[5], 
                                 noise = tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 2: Grow X stab CX: 5-2, 4-0

        self.circ.two_qubit_gate('CNOT', 
                            six_pack_col_list[4] + six_pack_col_list[5],
                            six_pack_col_list[0] + six_pack_col_list[2],
                            noise= tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] \
                            + six_pack_col_list[0] + six_pack_col_list[2]
        self._idle_depo(cleaness)
        self.circ.tick()


        # layer 3: Grow Z stab CX: 1-4, 3-5

        self.circ.two_qubit_gate('CNOT',
                                 six_pack_col_list[1] + six_pack_col_list[3],
                                 six_pack_col_list[4] + six_pack_col_list[5],
                                 noise= tq_noise)
        self.used_qubits = six_pack_col_list[1] + six_pack_col_list[3] \
                            + six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 4: Grow X and Z stab together CX: 0-1, 2-3
        self.circ.two_qubit_gate('CNOT',
                            six_pack_col_list[0] + six_pack_col_list[2],
                            six_pack_col_list[1] + six_pack_col_list[3],
                            noise= tq_noise)
        self.used_qubits = six_pack_col_list[0] + six_pack_col_list[2] \
                            + six_pack_col_list[1] + six_pack_col_list[3]
        self._idle_depo(cleaness)
        self.circ.tick()


        pass

    def color_to_rp5_no_growth(self, sc_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None


        six_packs = []
        for i in [-2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-1)])
        six_pack.append(self.qubit_network[(-1,1)])
        six_pack.append(self.qubit_network[(0,-1)])
        six_pack.append(self.qubit_network[(0,1)])
        six_pack.append(self.flag_network[(-1,-1)])
        six_pack.append(self.flag_network[(-1,1)])
        six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-2)])
        six_pack.append(self.qubit_network[(-1,2)])
        six_pack.append(self.qubit_network[(0,-2)])
        six_pack.append(self.qubit_network[(0,2)])
        six_pack.append(self.flag_network[(0,-1)])
        six_pack.append(self.flag_network[(0,1)])
        six_packs.append(six_pack)
        for i in [1,2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        
        six_pack_col_list = []
        for i in range(6):
            col_temp = []
            for six_pack in six_packs:
                col_temp.append(six_pack[i])
            six_pack_col_list.append(col_temp)

        # layer 4: Grow X and Z stab together CX: 0-1, 2-3
        self.circ.two_qubit_gate('CNOT',
                            six_pack_col_list[0] + six_pack_col_list[2],
                            six_pack_col_list[1] + six_pack_col_list[3],
                            noise=tq_noise)
        self.used_qubits = six_pack_col_list[0] + six_pack_col_list[2] \
                            + six_pack_col_list[1] + six_pack_col_list[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 3: Grow Z stab CX: 1-4, 3-5

        self.circ.two_qubit_gate('CNOT',
                                 six_pack_col_list[1] + six_pack_col_list[3],
                                 six_pack_col_list[4] + six_pack_col_list[5],
                                 noise= tq_noise)
        self.used_qubits = six_pack_col_list[1] + six_pack_col_list[3] \
                            + six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()

        self._idle_depo(cleaness)
        self.circ.tick()


        # layer 2: Grow X stab CX: 5-2, 4-0

        self.circ.two_qubit_gate('CNOT', 
                            six_pack_col_list[4] + six_pack_col_list[5],
                            six_pack_col_list[0] + six_pack_col_list[2],
                            noise= tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] \
                            + six_pack_col_list[0] + six_pack_col_list[2]
        self._idle_depo(cleaness)
        self.circ.tick()


        
        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]

        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list


        first_step = se_instruction.mid_cycle_steps[0]
        qubits_first_step = ctrl_dict[first_step] + targ_dict[first_step]


        # bell pair making

        self.circ.two_qubit_gate('CNOT', 
                                 six_pack_col_list[4],
                                 six_pack_col_list[5], 
                                 noise = tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] + qubits_first_step
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 1: Bell pair 
        self.circ.flag_measure_X(six_pack_col_list[4],
                                 self.stage,
                            noise=Z_error)
        self.circ.flag_measure(six_pack_col_list[5],
                               self.stage,
                            noise=X_error)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] + qubits_first_step
        self._pop_from_active(six_pack_col_list[4] + six_pack_col_list[5])
        self._idle_depo(cleaness)
        self.circ.tick()


        for six_pack in six_packs:
            self.circ.flag_detector(six_pack[4], self.stage, 
                        (six_pack[4].pos[0], six_pack[4].pos[1],
                         self.circ.time_temp, self.stage))
            self.circ.flag_detector(six_pack[5], self.stage, 
                        (six_pack[5].pos[0], six_pack[5].pos[1],
                         self.circ.time_temp, self.stage))
        
        self.stage_list.append(self.stage)
        self.stage += 1
        


    def color_to_rp5_new(self, cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        six_packs = []
        for i in [-2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-1)])
        six_pack.append(self.qubit_network[(-1,1)])
        six_pack.append(self.qubit_network[(0,-1)])
        six_pack.append(self.qubit_network[(0,1)])
        six_pack.append(self.flag_network[(-1,-1)])
        six_pack.append(self.flag_network[(-1,1)])
        six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-2)])
        six_pack.append(self.qubit_network[(-1,2)])
        six_pack.append(self.qubit_network[(0,-2)])
        six_pack.append(self.qubit_network[(0,2)])
        six_pack.append(self.flag_network[(0,-1)])
        six_pack.append(self.flag_network[(0,1)])
        six_packs.append(six_pack)
        for i in [1,2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        
        six_pack_col_list = []
        for i in range(6):
            col_temp = []
            for six_pack in six_packs:
                col_temp.append(six_pack[i])
            six_pack_col_list.append(col_temp)

        # layer 4: Grow X and Z stab together CX: 0-1, 2-3
        self.circ.two_qubit_gate('CNOT',
                            six_pack_col_list[0] + six_pack_col_list[2],
                            six_pack_col_list[1] + six_pack_col_list[3],
                            noise=tq_noise)
        self.used_qubits = six_pack_col_list[0] + six_pack_col_list[2] \
                            + six_pack_col_list[1] + six_pack_col_list[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 3: Grow Z stab CX: 1-4, 3-5

        self.circ.two_qubit_gate('CNOT',
                                 six_pack_col_list[1] + six_pack_col_list[3],
                                 six_pack_col_list[4] + six_pack_col_list[5],
                                 noise= tq_noise)
        self.used_qubits = six_pack_col_list[1] + six_pack_col_list[3] \
                            + six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 2: Grow X stab CX: 5-2, 4-0

        self.circ.two_qubit_gate('CNOT', 
                            six_pack_col_list[4] + six_pack_col_list[5],
                            six_pack_col_list[0] + six_pack_col_list[2],
                            noise= tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] \
                            + six_pack_col_list[0] + six_pack_col_list[2]
        self._idle_depo(cleaness)
        self.circ.tick()
        

        # bell pair making

        self.circ.two_qubit_gate('CNOT', 
                                 six_pack_col_list[4],
                                 six_pack_col_list[5], 
                                 noise = tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 1: Bell pair 
        self.circ.flag_measure_X(six_pack_col_list[4],
                                 self.stage,
                            noise=Z_error)
        self.circ.flag_measure(six_pack_col_list[5],
                               self.stage,
                            noise=X_error)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5]
        self._pop_from_active(six_pack_col_list[4] + six_pack_col_list[5])
        self._idle_depo(cleaness)
        self.circ.tick()


        for six_pack in six_packs:
            self.circ.flag_detector(six_pack[4], self.stage, 
                        (six_pack[4].pos[0], six_pack[4].pos[1],
                         self.circ.time_temp, self.stage))
            self.circ.flag_detector(six_pack[5], self.stage, 
                        (six_pack[5].pos[0], six_pack[5].pos[1],
                         self.circ.time_temp, self.stage))
        
        self.stage_list.append(self.stage)
        self.stage += 1
        

        pass

    
    def color_to_rp5_to_sc_bell(self, rp_name: str,
                                sc_name: str,
                                cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        six_packs = []
        for i in [-2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-1)])
        six_pack.append(self.qubit_network[(-1,1)])
        six_pack.append(self.qubit_network[(0,-1)])
        six_pack.append(self.qubit_network[(0,1)])
        six_pack.append(self.flag_network[(-1,-1)])
        six_pack.append(self.flag_network[(-1,1)])
        six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-2)])
        six_pack.append(self.qubit_network[(-1,2)])
        six_pack.append(self.qubit_network[(0,-2)])
        six_pack.append(self.qubit_network[(0,2)])
        six_pack.append(self.flag_network[(0,-1)])
        six_pack.append(self.flag_network[(0,1)])
        six_packs.append(six_pack)
        for i in [1,2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        
        six_pack_col_list = []
        for i in range(6):
            col_temp = []
            for six_pack in six_packs:
                col_temp.append(six_pack[i])
            six_pack_col_list.append(col_temp)

        # layer 4: Grow X and Z stab together CX: 0-1, 2-3
        self.circ.two_qubit_gate('CNOT',
                            six_pack_col_list[0] + six_pack_col_list[2],
                            six_pack_col_list[1] + six_pack_col_list[3],
                            noise=tq_noise)
        self.used_qubits = six_pack_col_list[0] + six_pack_col_list[2] \
                            + six_pack_col_list[1] + six_pack_col_list[3]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 3: Grow Z stab CX: 1-4, 3-5

        self.circ.two_qubit_gate('CNOT',
                                 six_pack_col_list[1] + six_pack_col_list[3],
                                 six_pack_col_list[4] + six_pack_col_list[5],
                                 noise= tq_noise)
        self.used_qubits = six_pack_col_list[1] + six_pack_col_list[3] \
                            + six_pack_col_list[4] + six_pack_col_list[5]
        self._idle_depo(cleaness)
        self.circ.tick()


        rp_code = self.code_name_dict[rp_name]
        sc_code = self.code_name_dict[sc_name]
        d_rp = rp_code.distance
        d_sc = sc_code.distance
        if d_rp >= d_sc:
            raise ValueError('rp code should be smaller')
        if d_sc - d_rp < 4:
            raise ValueError('not supported. growth too small')
        # step 1 reset 
        R_qubit_list = []
        RX_qubit_list = []
        bell_c_list = []
        bell_t_list = []
        y_1_larger = d_rp//2 +1
        for x in range(-y_1_larger, y_1_larger):
            bell_c_list.append(self.qubit_network[(x,-y_1_larger)])
            bell_t_list.append(self.qubit_network[(-x,y_1_larger)])
        for y in range(-y_1_larger+1,y_1_larger+1):
            bell_c_list.append(self.qubit_network[(-y_1_larger,y)])
            bell_t_list.append(self.qubit_network[(y_1_larger,-y)])
        for y in range(y_1_larger+1,d_sc//2+1):
            for x in range(-y,y):
                R_qubit_list.append(self.qubit_network[(x,-y)])
                R_qubit_list.append(self.qubit_network[(-x,y)])
        for x in range(y_1_larger+1,d_sc//2+1):
            for y in range(-x+1,x+1):
                RX_qubit_list.append(self.qubit_network[(-x,y)])
                RX_qubit_list.append(self.qubit_network[(x,-y)])

        code_temp = self.code_name_dict[sc_name]
        se_instruction = self.code_name_SE_instruction[sc_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])


        self.circ.reset(R_qubit_list + z_checks_in_circ + bell_t_list,
                        noise = ['X_ERROR', self.noise_probability])
        self.circ.reset_X(RX_qubit_list + x_checks_in_circ + bell_c_list,
                          noise = ['Z_ERROR', self.noise_probability])
        self._add_to_active(R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ
                            +bell_t_list+bell_c_list)
        self.used_qubits = R_qubit_list+z_checks_in_circ+RX_qubit_list+x_checks_in_circ \
                            +bell_t_list+bell_c_list
        self._idle_depo()
        self.circ.tick()


        # layer 2: Grow X stab CX: 5-2, 4-0

        self.circ.two_qubit_gate('CNOT', 
                            six_pack_col_list[4] + six_pack_col_list[5] + bell_c_list,
                            six_pack_col_list[0] + six_pack_col_list[2] + bell_t_list,
                            noise= tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] \
                            + six_pack_col_list[0] + six_pack_col_list[2] \
                            + bell_t_list + bell_c_list
        self._idle_depo(cleaness)
        self.circ.tick()



        ctrl_dict = {}
        targ_dict = {}
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass
            ctrl_dict[i] = ctrl_list
            targ_dict[i] = targ_list


        first_step = se_instruction.mid_cycle_steps[0]

        

        # bell pair making

        self.circ.two_qubit_gate('CNOT', 
                                 six_pack_col_list[4] + ctrl_dict[first_step],
                                 six_pack_col_list[5] + targ_dict[first_step], 
                                 noise = tq_noise)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5] \
                            + ctrl_dict[first_step] + targ_dict[first_step]
        self._idle_depo(cleaness)
        self.circ.tick()

        # layer 1: Bell pair 
        self.circ.flag_measure_X(six_pack_col_list[4],
                                 self.stage,
                            noise=Z_error)
        self.circ.flag_measure(six_pack_col_list[5],
                               self.stage,
                            noise=X_error)
        self.used_qubits = six_pack_col_list[4] + six_pack_col_list[5]
        self._pop_from_active(six_pack_col_list[4] + six_pack_col_list[5])
        self._idle_depo(cleaness)
        self.circ.tick()


        for six_pack in six_packs:
            self.circ.flag_detector(six_pack[4], self.stage, 
                        (six_pack[4].pos[0], six_pack[4].pos[1],
                         self.circ.time_temp, self.stage))
            self.circ.flag_detector(six_pack[5], self.stage, 
                        (six_pack[5].pos[0], six_pack[5].pos[1],
                         self.circ.time_temp, self.stage))
        
        self.stage_list.append(self.stage)
        self.stage += 1
        

        for i in se_instruction.mid_cycle_steps:
            if i != first_step:
                self.circ.two_qubit_gate('CX', ctrl_dict[i],
                                            targ_dict[i],
                                            noise=tq_noise)
                self.used_qubits = ctrl_dict[i]+targ_dict[i]
                self._idle_depo()
                self.circ.tick()

        self.circ.measure(z_checks_in_circ,noise=X_error)
        self.circ.measure_X(x_checks_in_circ,noise=Z_error)
        self.used_qubits = z_checks_in_circ + x_checks_in_circ
        self._idle_depo()
        self._pop_from_active(z_checks_in_circ+x_checks_in_circ)
        self.circ.tick()

        # anc positions in rp code:
        anc_pos_set = set()
        for check in rp_code.x_check_collection + rp_code.z_check_collection:
            anc_pos_set.add(check.pos)
        
        # detector specification:
        for check in sc_code.x_check_collection + sc_code.z_check_collection:
            qubit = self.qubit_network[check.pos]
            if np.abs(qubit.pos[0]) < d_rp // 2 and \
               np.abs(qubit.pos[1]) < d_rp // 2:
                self.circ.detector_qt_list([qubit, qubit],
                                [self.circ.time_temp-1, 
                                 self.circ.time_temp],
                         (qubit.pos[0],qubit.pos[1],self.circ.time_temp))
            elif qubit.pos in anc_pos_set:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit,qubit_inv],
                                           [self.circ.time_temp-1,
                                            self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) < (d_rp // 2) + 1 and \
                 qubit.pos[1] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) < (d_rp // 2) + 1 and \
                 qubit.pos[0] == -(d_rp//2) - 1.5 and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile:
                qubit_inv_pos = (-qubit.pos[0],-qubit.pos[1])
                qubit_inv = self.qubit_network[qubit_inv_pos]
                self.circ.detector_qt_list([qubit,qubit_inv],
                                           [self.circ.time_temp,
                                            self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[1]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[1]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[1]) > np.abs(qubit.pos[0]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 1: # blue tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))
            elif np.abs(qubit.pos[0]) > d_rp//2 + 2 and \
                 np.abs(qubit.pos[0]) < d_sc//2 + 1 and \
                 np.abs(qubit.pos[0]) > np.abs(qubit.pos[1]) and \
                 (qubit.pos[0]-qubit.pos[1]) % 2 == 0: # red tile
                self.circ.detector_qt_list([qubit],
                                           [self.circ.time_temp],
                     (qubit.pos[0], qubit.pos[1], self.circ.time_temp))            
        self.circ.clock_plus_one()
        self.circ.tick()        

        pass



    def d_5_color_T_check(self, rp5_name: str, 
                          cleaness: bool = False) -> None:
        if cleaness == False:
            sq_noise = ['DEPOLARIZE1',self.noise_probability]
            tq_noise = ['DEPOLARIZE2',self.noise_probability]
            X_error = ['X_ERROR',self.noise_probability]
            Z_error = ['Z_ERROR',self.noise_probability]
        else:
            sq_noise = None
            tq_noise = None
            X_error = None
            Z_error = None

        code = self.code_name_dict[rp5_name]

        six_packs = []

        for i in [-2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-1)])
        six_pack.append(self.qubit_network[(-1,1)])
        six_pack.append(self.qubit_network[(0,-1)])
        six_pack.append(self.qubit_network[(0,1)])
        six_pack.append(self.flag_network[(-1,-1)])
        six_pack.append(self.flag_network[(-1,1)])
        six_packs.append(six_pack)
        six_pack = []
        six_pack.append(self.qubit_network[(-1,-2)])
        six_pack.append(self.qubit_network[(-1,2)])
        six_pack.append(self.qubit_network[(0,-2)])
        six_pack.append(self.qubit_network[(0,2)])
        six_pack.append(self.flag_network[(0,-1)])
        six_pack.append(self.flag_network[(0,1)])
        six_packs.append(six_pack)
        for i in [1,2]:
            six_pack = []
            six_pack.append(self.qubit_network[(i,-2)])
            six_pack.append(self.qubit_network[(i,2)])
            six_pack.append(self.qubit_network[(i,-1)])
            six_pack.append(self.qubit_network[(i,1)])
            six_pack.append(self.flag_network[(i,-1)])
            six_pack.append(self.flag_network[(i,1)])
            six_packs.append(six_pack)

        shadow_qubits: list[Qubit] = []
        real_qubits = []
        image_qubits = []
        for i in [-2,-1,0,1,2]:
            if i == -1:
                shadow_qubits.append(self.qubit_network[(-1.5,0.5)])
                image_qubits.append(self.qubit_network[(-1,1)])
                shadow_qubits.append(self.qubit_network[(-0.5,0.5)])
                image_qubits.append(self.qubit_network[(0,1)])
            elif i == 0:
                shadow_qubits.append(self.qubit_network[(-1.5,1.5)])
                image_qubits.append(self.qubit_network[(-1,2)])
                shadow_qubits.append(self.qubit_network[(-0.5,1.5)])
                image_qubits.append(self.qubit_network[(0,2)])
            else:
                shadow_qubits.append(self.qubit_network[(i-0.5,0.5)])
                image_qubits.append(self.qubit_network[(i,1)])
                shadow_qubits.append(self.qubit_network[(i-0.5,1.5)])
                image_qubits.append(self.qubit_network[(i,2)])
        for i in [-2,-1,0,1,2]:
            shadow_qubits.append(self.qubit_network[(i-0.5,-0.5)])
            image_qubits.append(self.qubit_network[(i,0)])
        for six_pack in six_packs:
            for qubit in six_pack:
                real_qubits.append(qubit)
        for i in [-2,-1,0,1,2]:
            real_qubits.append(self.qubit_network[(i,0)])
        

        self.circ.single_qubit_gate('S_DAG', real_qubits, sq_noise)
        self.used_qubits = real_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        self.circ.reset_X(shadow_qubits, Z_error)
        self._add_to_active(shadow_qubits)
        self.used_qubits = shadow_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        

        self.circ.two_qubit_gate('CNOT', shadow_qubits,
                                 image_qubits, tq_noise)
        self.used_qubits = shadow_qubits + image_qubits
        self._idle_depo(cleaness)
        self.circ.tick()
        

        ctrl_qubits_dict = {}
        targ_qubits_dict = {}

        ctrl_qubits_dict[4] = [self.qubit_network[(i,0)] 
                                for i in range(-2,3)]
        targ_qubits_dict[4] = [six_packs[i][0] 
                                for i in range(0,5)]
        ctrl_qubits_dict[5] = [self.qubit_network[(i,0)]
                                for i in range(-2,3)] \
                            + [six_packs[i][0]
                                for i in range(0,5)]
        targ_qubits_dict[5] = [six_packs[i][3]
                                for i in range(0,5)] \
                            + [six_packs[i][1]
                                for i in range(0,5)]
        ctrl_qubits_dict[6] = [self.qubit_network[(i,0)]
                                for i in range(-2,3)] \
                            + [six_packs[i][1]
                                for i in range(0,5)] \
                            + [six_packs[i][3]
                                for i in range(0,5)]
        targ_qubits_dict[6] = [six_packs[i][2]
                                for i in range(0,5)] \
                            + [six_packs[i][4]
                                for i in range(0,5)] \
                            + [six_packs[i][5]
                                for i in range(0,5)]
        ctrl_qubits_dict[1] = [self.qubit_network[(0,0)]]
        targ_qubits_dict[1] = [self.qubit_network[(-1,0)]]
        ctrl_qubits_dict[2] = [self.qubit_network[(0,0)],
                               self.qubit_network[(-1,0)]]
        targ_qubits_dict[2] = [self.qubit_network[(1,0)],
                               self.qubit_network[(-2,0)]]
        ctrl_qubits_dict[3] = [self.qubit_network[(0,0)]]
        targ_qubits_dict[3] = [self.qubit_network[(2,0)]]

        for t in [6,5,4,3,2,1]:
            self.circ.two_qubit_gate('CNOT', ctrl_qubits_dict[t],
                                    targ_qubits_dict[t], tq_noise)
            self.used_qubits = ctrl_qubits_dict[t] + targ_qubits_dict[t]
            self._idle_depo(cleaness)
            self.circ.tick()
        
        self.circ.flag_measure_X([self.qubit_network[(0,0)]],
                                 self.stage,
                            Z_error)
        self.used_qubits = [self.qubit_network[(0,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()

        self.circ.flag_detector(self.qubit_network[(0,0)],
                                   self.stage,
                                    (0,0,self.circ.time_temp,self.stage))
        self.circ.reset_X([self.qubit_network[(0,0)]], Z_error)
        self.used_qubits = [self.qubit_network[(0,0)]]
        self._idle_depo(cleaness)
        self.circ.tick()
        self.stage_list.append(self.stage)
        self.stage += 1


        for t in [1,2,3,4,5,6]:
            self.circ.two_qubit_gate('CNOT', ctrl_qubits_dict[t],
                                    targ_qubits_dict[t], tq_noise)
            self.used_qubits = ctrl_qubits_dict[t] + targ_qubits_dict[t]
            self._idle_depo(cleaness)
            self.circ.tick()

        
        self.circ.two_qubit_gate('CNOT', shadow_qubits,
                                 image_qubits, tq_noise)
        self.used_qubits = shadow_qubits + image_qubits
        self._idle_depo(cleaness)
        self.circ.tick()
        

        self.circ.flag_measure_X(shadow_qubits, self.stage, Z_error)
        self.used_qubits = shadow_qubits 
        self._idle_depo(cleaness)
        self._pop_from_active(shadow_qubits)
        self.circ.tick()

        self.circ.single_qubit_gate('S', real_qubits, sq_noise)
        self.used_qubits = real_qubits
        self._idle_depo(cleaness)
        self.circ.tick()

        for shadow in shadow_qubits:
            self.circ.flag_detector(shadow,
                                    self.stage,
                (shadow.pos[0], shadow.pos[1], self.circ.time_temp, self.stage))
        
        self.stage_list.append(self.stage)
        self.stage += 1

        pass



        

class SE_instruction_set():

    def __init__(self, code_name: str,
                 mid_cycle_steps: list,
                 z_ctrl_dict: dict[int,list[Qubit]],
                 z_targ_dict: dict[int,list[Qubit]],
                 x_ctrl_dict: dict[int,list[Qubit]],
                 x_targ_dict: dict[int,list[Qubit]]
                 ):
        
        self.code_name = code_name
        self.mid_cycle_steps = mid_cycle_steps
        self.z_ctrl_dict = z_ctrl_dict
        self.z_targ_dict = z_targ_dict
        self.x_ctrl_dict = x_ctrl_dict
        self.x_targ_dict = x_targ_dict




def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors


