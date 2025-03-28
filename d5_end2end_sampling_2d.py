from t_gadget import *
from so_sampler_2d import *
from csv_so_processor import *


def main():
    d_rp3 = 3
    d_rp5 = 5
    d_sc = 11
    rp3_name = 'rp' + str(d_rp3)
    rp5_name = 'rp' + str(d_rp5)
    sc_name = 'sc' + str(d_sc)
    rp3_code = RP_surface_code_fresh(d_rp3, rp3_name)
    rp5_code = RP_surface_code_fresh(d_rp5, rp5_name)
    sc_code = Rotated_surface_code(d_sc, sc_name)
    n_rounds = d_rp5 + 1
    path_name = './circuit_garage/'
    circ_name =  path_name + 'rp_' + str(d_rp3) \
                        + '_rp_' + str(d_rp5) \
                        + '_sc_' + str(d_sc) + '_end2end_' \
                        + str(n_rounds) + '_full_rds'  + '.stim'
    circ_svg_name = './circuit_gallery/' + 'rp_' + str(d_rp3) \
                        + '_rp_' + str(d_rp5) \
                        + '_sc_' + str(d_sc) + '_end2end_' \
                        + str(n_rounds) + '_full_rds'  + '.svg'
    result_path = './sample_results/' + 'rp_' + str(d_rp3) \
                    + '_rp_' + str(d_rp5) \
                    + '_sc_' + str(d_sc) + '_end2end_' \
                    + str(n_rounds) + '_full_rds' + '_so2d.csv'
    
    result_temp_path = './sample_results/' + 'rp_' + str(d_rp3) \
                    + '_rp_' + str(d_rp5) \
                    + '_sc_' + str(d_sc) + '_end2end_' \
                    + str(n_rounds) + '_full_rds' + '_so2d_temp.csv'
    result_name_no_append = 'rp_' + str(d_rp3) \
                    + '_rp_' + str(d_rp5) \
                    + '_sc_' + str(d_sc) + '_end2end_' \
                    + str(n_rounds) + '_full_rds' + '_so2d'

    sc_code.logic_x_selection(0)
    sc_code.logic_z_selection(0)


    noise = 0.001

    cleaness = False
    detector_val = False

    ft_circuit = Meta_Circuit([rp3_code,rp5_code,sc_code],circ_name,noise)

    ft_circuit.Y_to_rp3_sign_correction(rp3_name, cleaness=cleaness,
                                        detector_val=detector_val)

    for i in range(1):
        ft_circuit.SE_round_z_corrected_rp3(rp3_name, cleaness=cleaness, 
                                        post=True, z_sign_corrected=True,
                                        detector_val=detector_val)




    ft_circuit.rp3_to_color_compact(rp3_name,cleaness=cleaness)
    ft_circuit.d3_color_T_check_compact(rp3_name, cleaness=cleaness,
                                        detector_val=detector_val)
    ft_circuit.color_to_rp3_to_rp5(rp3_name, rp5_name)
    ft_circuit.SE_round_z_corrected(rp5_name, cleaness=cleaness, post=True,
                                    z_sign_corrected=True)
    ft_circuit.rp5_to_color_new()
    ft_circuit.d_5_color_T_check(rp5_name)
    ft_circuit.color_to_rp5_to_sc_bell(rp5_name,sc_name)


    for i in range(d_rp5):
        ft_circuit.SE_round(sc_name)

    ft_circuit.magic_Y_measurement(sc_name)

    # ft_circuit.check_circuit_distance()


    stim_circ = ft_circuit.get_stim_circuit()



    with open(circ_svg_name,'w') as file:
        out = stim_circ.diagram('detslice-with-ops-svg')
        file.write(out.__str__())


    max_shots = 2_000_000_000_000
    iterations = 4_000
    shots_per_iter = max_shots // iterations


    sinter_task = sinter.Task(circuit=stim_circ)

    res_temp = SO_2d(file_name_no_append=result_name_no_append)
    res_temp.clear_cache()

    for iter in range(iterations):
        print(iter,'/',iterations)
        samples = sinter.collect(
            num_workers=50,
            max_shots=shots_per_iter,
            # max_errors=1000,
            tasks=[sinter_task],
            decoders=['SO5'],
            custom_decoders={'SO5':SO5Sampler2D()},
            save_resume_filepath=result_temp_path,
            print_progress=False
        )
        res_temp = SO_2d(file_name_no_append=result_name_no_append)
        res_temp.update()



# NOTE: This is actually necessary! If the code inside 'main()' was at the
# module level, the multiprocessing children spawned by sinter.collect would
# also attempt to run that code.
if __name__ == '__main__':
    main()


