import stim
import sinter
import numpy as np
from typing import Optional

class Stage():

    def __init__(self, 
                 stage_name: str,
                 start: list[str], 
                 end: list[str]):
        '''
        If both start and end are non-empty, 
        '''
        self.stage_name = stage_name
        self.start = start
        self.end = end
        self.start_opline: Optional[int] = None
        self.end_opline: Optional[int] = None
        self.active_qubits: Optional[int] = 0
        self.det_stage: Optional[int] = None
        self.survival_rate: Optional[float] = None

class Circ_parsor():

    def __init__(self, circ_file_path: str):
        self.circ = stim.Circuit.from_file(circ_file_path)

        self.break_points: list[tuple[int,str]] = []
        self.stages: list[Stage] = []
        self.post_stages: list[Stage] = []
        self.unpost_stages: list[Stage] = []
        self.gap_survival_rate: Optional[float] = None
        self.before_final_sr: Optional[float] = None
        pass

    def config_stages(self) -> None:
        pass

    def load_from_discard_tests(self, result_file_path: str) -> None:
        stats = sinter.read_stats_from_csv_files(result_file_path)
        if len(stats) != 1:
            raise NotImplementedError('multiple task results in a single csv file')
        stats = stats[0]
        shots = stats.shots
        discard_temp = 0
        det_stages = sorted([int(key[1:]) for key in stats.custom_counts])
        survival_rates = []
        for det_stage in det_stages:
            discard_temp += stats.custom_counts['D'+str(det_stage)]
            survival_rates.append((shots-discard_temp)/shots)

        self.before_final_sr = survival_rates[-2]
        print(discard_temp,stats.discards)
        assert discard_temp == stats.discards
        for stage in self.post_stages:
            if stage.det_stage == 0:
                stage.survival_rate = 1
            else:
                stage.survival_rate = survival_rates[stage.det_stage-1]
        
        for stage in self.stages:
            if stage not in self.post_stages:
                stage.survival_rate = self.before_final_sr
        # for stage in self.unpost_stages:
        #     stage.survival_rate = self.post_stages[-1].
        pass

    def active_qubits_calc(self) -> None:
        trivial_ops = {'DEPOLARIZE1','DEPOLARIZE2','X_ERROR','Z_ERROR','I',
                       'QUBIT_COORDS','TICK','DETECTOR','OBSERVABLE_INCLUDE'}
        RM_ops = {'R','RX','M','MX','MPP'}
        current_stage = self.stages[0]
        current_stage_bd = current_stage.start + current_stage.end
        current_start = [i for i in current_stage.start]
        current_end = [i for i in current_stage.end]
        active_qubit_ids = set()
        current_stage_index = 0
        for instruction in self.circ:
            if instruction.name not in RM_ops:
                continue
            try:
                current_stage_bd.remove(instruction.name)
                #print(f"stage {current_stage.stage_name} removed {instruction.name}")
            except ValueError:
                print(f"stage {current_stage.stage_name} did not expect {instruction.name}")
                raise
            gate_targets = instruction.targets_copy()
            if instruction.name in ['R','RX']:
                for gate_target in gate_targets:
                    active_qubit_ids.add(gate_target.value)
                current_start.remove(instruction.name)
                if len(current_start) == 0:
                    current_stage.active_qubits = len(active_qubit_ids)
            elif instruction.name in ['MX','M']:
                for gate_target in gate_targets:
                    active_qubit_ids.remove(gate_target.value)
            elif instruction.name == 'MPP':
                break
            if len(current_stage_bd) == 0:
                current_stage_index += 1
                current_stage = self.stages[current_stage_index]
                current_stage_bd = current_stage.start + current_stage.end
                current_start = [i for i in current_stage.start]
                if len(current_stage.start) == 0:
                    current_stage.active_qubits = len(active_qubit_ids)
        

        pass


    def test_discard(self, result_fp1, result_fp2) -> None:
        pass

    def volume_calc(self, final_survival_rate: float) -> int:
        vol = 0
        for stage in self.stages:
            if stage in self.post_stages:
                if 'Morph' in stage.stage_name:
                    if 'Back&Grow' in stage.stage_name:
                        vol += stage.active_qubits * stage.survival_rate*1.5
                    else:
                        vol += stage.active_qubits * stage.survival_rate * 0.5
                elif 'PerfMeas' in stage.stage_name:
                    pass
                else:
                    vol += stage.active_qubits * stage.survival_rate
            else:
                if 'PerfMeas' in stage.stage_name:
                    vol += self.stages[-2].active_qubits * final_survival_rate
                else:
                    vol += stage.active_qubits * self.before_final_sr
        return round(vol/final_survival_rate)


class Circ_rp3_T_ungrown(Circ_parsor):

    def __init__(self, circ_file_path: str):
        super().__init__(circ_file_path)
    
    def config_stages(self) -> None:
        det_stage = 0
        # injection :
        injection_stage = Stage('Injection',
                                ['R','RX','RX','R','RX'],['M','MX'])
        injection_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(injection_stage)
        self.post_stages.append(injection_stage)
        # stab meas rp3
        se_stage = Stage('SE',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp3
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX'])
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back to rp3
        morphback_stage = Stage('MorphBack',[],['M','MX'])
        morphback_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(morphback_stage)
        self.post_stages.append(morphback_stage)
        # magic measurement
        magic_meas_stage = Stage('PerfMeas',[],['MPP'])
        magic_meas_stage.det_stage = det_stage
        self.stages.append(magic_meas_stage)
        self.post_stages.append(magic_meas_stage)

        for i in reversed(range(len(self.post_stages))):
            stage = self.post_stages[i]
            if stage.det_stage is None:
                for j in range(i+1,len(self.post_stages)):
                    stage_j = self.post_stages[j]
                    if stage_j.det_stage is not None:
                        stage.det_stage = stage_j.det_stage
                        break
                



        pass



class Circ_rp3_T_end2end(Circ_parsor):

    def __init__(self, circ_file_path: str, 
                 sc_d: int, n_rounds: int):
        super().__init__(circ_file_path)
        self.sc_d = sc_d
        self.n_rounds = n_rounds
    
    def config_stages(self) -> None:
        det_stage = 0
        # injection :
        injection_stage = Stage('Injection',
                                ['R','RX','RX','R','RX'],['M','MX'])
        injection_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(injection_stage)
        self.post_stages.append(injection_stage)
        # stab meas rp3
        se_stage = Stage('SE',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp3
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX']) 
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back and grow
        morphback_and_grow_stage = Stage('MorphBack&Grow',['R','RX'],['M',
                                                    'MX','M','MX'])
        morphback_and_grow_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(morphback_and_grow_stage)
        self.post_stages.append(morphback_and_grow_stage)
        # SE rounds: 
        for i in range(self.n_rounds):
            se_stage = Stage('SE_SC-'+str(i+1),['R','RX'],
                             ['M','MX'])
            self.stages.append(se_stage)



        # magic measurement
        magic_meas_stage = Stage('PerfMeas',[],['MPP'])
        self.stages.append(magic_meas_stage)

        for i in reversed(range(len(self.post_stages))):
            stage = self.post_stages[i]
            if stage.det_stage is None:
                for j in range(i+1,len(self.post_stages)):
                    stage_j = self.post_stages[j]
                    if stage_j.det_stage is not None:
                        stage.det_stage = stage_j.det_stage
                        break
                



        pass


class Circ_rp5_T_ungrown(Circ_parsor):

    def __init__(self, circ_file_path):
        super().__init__(circ_file_path)
    
    def config_stages(self) -> None:
        
        det_stage = 0
        # injection :
        injection_stage = Stage('Injection',
                                ['R','RX','RX','R','RX'],['M','MX'])
        injection_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(injection_stage)
        self.post_stages.append(injection_stage)
        # stab meas rp3
        se_stage = Stage('SE',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp3
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX'])
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back to rp3 and grow to rp5
        morphback_stage = Stage('MorphBack&Grow',['R','RX'],
                                ['M','MX','M','MX'])
        morphback_stage.det_stage = det_stage + 1
        det_stage += 2
        self.stages.append(morphback_stage)
        self.post_stages.append(morphback_stage)
        # SE-RP5
        se_stage = Stage('SE_RP5',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp5
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX'])
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back to rp5
        morphback_stage = Stage('MorphBack',[],['M','MX'])
        morphback_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(morphback_stage)
        self.post_stages.append(morphback_stage)
        # magic measurement
        magic_meas_stage = Stage('PerfMeas',[],['MPP'])
        magic_meas_stage.det_stage = det_stage
        self.stages.append(magic_meas_stage)
        self.post_stages.append(magic_meas_stage)

        for i in reversed(range(len(self.post_stages))):
            stage = self.post_stages[i]
            if stage.det_stage is None:
                for j in range(i+1,len(self.post_stages)):
                    stage_j = self.post_stages[j]
                    if stage_j.det_stage is not None:
                        stage.det_stage = stage_j.det_stage
                        break
                







class Circ_rp5_T_end2end(Circ_parsor):

    def __init__(self, circ_file_path,
                 sc_d: int, n_rounds: int):
        super().__init__(circ_file_path)
        self.sc_d = sc_d
        self.n_rounds = n_rounds
    
    def config_stages(self) -> None:
        
        det_stage = 0
        # injection :
        injection_stage = Stage('Injection',
                                ['R','RX','RX','R','RX'],['M','MX'])
        injection_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(injection_stage)
        self.post_stages.append(injection_stage)
        # stab meas rp3
        se_stage = Stage('SE',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp3
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX'])
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back to rp3 and grow to rp5
        morphback_stage = Stage('MorphBack&Grow',['R','RX'],
                                ['M','MX','M','MX'])
        morphback_stage.det_stage = det_stage + 1
        det_stage += 2
        self.stages.append(morphback_stage)
        self.post_stages.append(morphback_stage)
        # SE-RP5
        se_stage = Stage('SE_RP5',['R','RX'],['M','MX'])
        se_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(se_stage)
        self.post_stages.append(se_stage)
        # morph to srp5
        morphto_stage = Stage('MorphTo',['R','RX'],[])
        self.stages.append(morphto_stage)
        self.post_stages.append(morphto_stage)
        # T check
        t_check_stage = Stage('T_check',['RX'],['MX'])
        t_check_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_stage)
        self.post_stages.append(t_check_stage)
        # T check reverse
        t_check_inv_stage = Stage('T_check',['RX'],['MX'])
        t_check_inv_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(t_check_inv_stage)
        self.post_stages.append(t_check_inv_stage)
        # morph back and grow
        morphback_and_grow_stage = Stage('MorphBack&Grow',['R','RX'],['M',
                                                    'MX','M','MX'])
        morphback_and_grow_stage.det_stage = det_stage
        det_stage += 1
        self.stages.append(morphback_and_grow_stage)
        self.post_stages.append(morphback_and_grow_stage)

        # SE rounds: 
        for i in range(self.n_rounds):
            se_stage = Stage('SE_SC-'+str(i+1),['R','RX'],
                             ['M','MX'])
            self.stages.append(se_stage)

        # magic measurement
        magic_meas_stage = Stage('PerfMeas',[],['MPP'])
        self.stages.append(magic_meas_stage)

        for i in reversed(range(len(self.post_stages))):
            stage = self.post_stages[i]
            if stage.det_stage is None:
                for j in range(i+1,len(self.post_stages)):
                    stage_j = self.post_stages[j]
                    if stage_j.det_stage is not None:
                        stage.det_stage = stage_j.det_stage
                        break
                

