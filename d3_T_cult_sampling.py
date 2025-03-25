from t_gadget import *
from perfectionist_sampler import *
from csv_so_processor import *


def main():
    d_rp3 = 3
    d_sc = 7
    rp3_name = 'rp' + str(d_rp3)
    sc_name = 'sc' + str(d_sc)
    rp3_code = RP_surface_code_fresh(d_rp3, rp3_name)
    sc_code = Rotated_surface_code(d_sc, sc_name)
    path_name = './circuit_garage/'
    circ_name =  path_name + 'rp_' + str(d_rp3) + '_T_cult.stim'
    circ_svg_name = './circuit_gallery/' + 'rp_' + str(d_rp3) + '_T_cult.svg' 
    result_temp_path = './sample_results/' + 'rp_' + str(d_rp3) + '_T_cult_temp.csv'
    result_name_no_append = 'rp_' + str(d_rp3) + '_T_cult'
    # sc_code.qubit_displacement(-((d_sc-d_rp3)//2),
    #                            -((d_sc-d_rp3)//2))
    sc_code.logic_x_selection(0)
    sc_code.logic_z_selection(0)


    noise = 0.001

    cleaness = False
    detector_val = False

    ft_circuit = Meta_Circuit([rp3_code,sc_code],circ_name,noise)

    ft_circuit.Y_to_rp3_sign_correction(rp3_name, cleaness=cleaness,
                                        detector_val=detector_val)

    for i in range(1):
        ft_circuit.SE_round_z_corrected_rp3(rp3_name, cleaness=cleaness, 
                                        post=True, z_sign_corrected=True,
                                        detector_val=detector_val)




    ft_circuit.rp3_to_color_compact(rp3_name,cleaness=cleaness)
    ft_circuit.d3_color_T_check_compact(rp3_name, cleaness=cleaness,
                                        detector_val=detector_val)
    ft_circuit.color_to_rp3_compact_no_growth(rp3_name,sc_name)



    ft_circuit.magic_Y_measurement_post_selection(rp3_name)



    stim_circ = ft_circuit.get_stim_circuit()




    with open(circ_svg_name,'w') as file:
        out = stim_circ.diagram('detslice-with-ops-svg')
        file.write(out.__str__())



    max_shots = 2_000_000_000
    iterations = 5
    shots_per_iter = max_shots // iterations


    sinter_task = sinter.Task(circuit=stim_circ)
    res_temp = Perf(file_name_no_append=result_name_no_append)
    res_temp.clear_cache()

    for _ in range(iterations):
        samples = sinter.collect(
            num_workers=50,
            max_shots=shots_per_iter,
            # max_errors=1000,
            tasks=[sinter_task],
            decoders=['perfect'],
            custom_decoders={'perfect':PerfectionistSampler()},
            save_resume_filepath=result_temp_path,
            print_progress=True
        )
        res_temp = Perf(file_name_no_append=result_name_no_append)
        res_temp.update()

        



# NOTE: This is actually necessary! If the code inside 'main()' was at the
# module level, the multiprocessing children spawned by sinter.collect would
# also attempt to run that code.
if __name__ == '__main__':
    main()


