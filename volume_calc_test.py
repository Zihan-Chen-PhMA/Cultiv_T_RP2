import stim
from volume_calc import *

filename = 'rp_3_T_cult'

circ_path = './circuit_garage/' + filename + '.stim'
res_path = './sample_results/' + filename + '_combined.csv'
discard_res_path = './sample_results/' + filename + '_combined.csv'


vol_helper = Circ_rp3_T_ungrown(circ_path)
vol_helper.config_stages()
vol_helper.active_qubits_calc()
vol_helper.load_from_discard_tests(res_path)

for stage in vol_helper.stages:
    print(stage.stage_name,stage.active_qubits,stage.survival_rate,stage.det_stage)


filename = 'rp_3_sc_7_end2end_4_full_rds'
circ_path = './circuit_garage/' + filename + '.stim'
res_path = './sample_results/' + filename + '_so2d_combined.csv'
sc_d = 7
n_rounds = 2

vol_helper = Circ_rp3_T_end2end(circ_path,sc_d,n_rounds)
vol_helper.config_stages()
vol_helper.active_qubits_calc()
vol_helper.load_from_discard_tests(discard_res_path)

for stage in vol_helper.stages:
    print(stage.stage_name, stage.active_qubits, stage.survival_rate, stage.det_stage)

print(vol_helper.volume_calc(1/2.3),vol_helper.volume_calc(1/2.5))





filename = 'rp_3_rp_5_T_cult'
circ_path = './circuit_garage/' + filename + '.stim'
# res_path = './sample_results/' + filename + '_so2d_combined.csv'
discard_res_path = './sample_results/' + filename + '_combined.csv'
sc_d = 7
n_rounds = 2

vol_helper = Circ_rp5_T_ungrown(circ_path)
vol_helper.config_stages()
vol_helper.active_qubits_calc()
vol_helper.load_from_discard_tests(discard_res_path)

for stage in vol_helper.stages:
    print(stage.stage_name, stage.active_qubits, stage.survival_rate, stage.det_stage)

print(vol_helper.volume_calc(1/10))



filename = 'rp_3_rp_5_sc_11_end2end_6_full_rds'
circ_path = './circuit_garage/' + filename + '.stim'
res_path = './sample_results/' + filename + '_so2d_combined.csv'
sc_d = 11
n_rounds = 5

vol_helper = Circ_rp5_T_end2end(circ_path,sc_d,n_rounds)
vol_helper.config_stages()
vol_helper.active_qubits_calc()
vol_helper.load_from_discard_tests(discard_res_path)

for stage in vol_helper.stages:
    print(stage.stage_name, stage.active_qubits, stage.survival_rate, stage.det_stage)

print(vol_helper.volume_calc(1/15))
