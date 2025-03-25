import sinter
import numpy as np
import os.path as path


class RES_Handler():

    def __init__(self, 
                 file_name_no_append: str,
                 pathdir: str):
        self.filename_temp = pathdir + file_name_no_append + '_temp.csv'
        self.filename_save = pathdir + file_name_no_append + '_combined.csv'
        # if (path.isfile(self.filename_temp) or path.isfile(self.filename_save)) == False:
        #     raise ValueError('No results present.')
        self.shots = 0
        self.discards = 0
        self.errors = 0
        pass

    def clear_cache(self) -> None:
        if path.isfile(self.filename_temp):
            with open(self.filename_temp,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')

    def update(self) -> None:
        if path.isfile(self.filename_temp) == False:
            raise ValueError('no new data present')
        stats_temp = sinter.read_stats_from_csv_files(self.filename_temp)
        if len(stats_temp) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        stats_temp = stats_temp[0]
        if path.isfile(self.filename_save) == False:
            with open(self.filename_save,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')
                file.write(stats_temp.to_csv_line())
                file.write('\n')
        else:
            existing_stats = sinter.read_stats_from_csv_files(self.filename_save)
            if len(existing_stats) != 1:
                raise ValueError('More than one task result is stored. Unsupported behavior for now.')
            existing_stats = existing_stats[0]
            updated_stats = stats_temp + existing_stats
            with open(self.filename_save,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')
                file.write(updated_stats.to_csv_line())
                file.write('\n')
        
        with open(self.filename_temp,'w') as file:
            file.write(sinter.CSV_HEADER)
            file.write('\n')

    def rate(self, hits_list: list[int], shots_list: list[int]) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=1000).best)
        return ret_list

    def rate_ceil_binfit(self, hits_list: list[int],
                         shots_list: list[int],
                         param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).high)
        return ret_list

    def rate_floor_binfit(self, hits_list: list[int],
                          shots_list: list[int],
                          param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).low)
        return ret_list
    
    def read_through_custom_counts(self) -> None:
        if path.isfile(self.filename_save) == False:
            raise ValueError('no result found')
        self.stats = sinter.read_stats_from_csv_files(self.filename_save)
        if len(self.stats) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        self.stats = self.stats[0]
        self.shots = self.stats.shots
        self.discards = self.stats.discards
        self.errors = self.stats.errors
    

class Perf(RES_Handler):

    def __init__(self,
                 file_name_no_append: str,
                 pathdir: str = './sample_results/'):
        super().__init__(file_name_no_append, pathdir)
        


class Gap(RES_Handler):

    def __init__(self, pathdir: str = './msc_paper_results/'):
        super().__init__('', pathdir)
        self.filename_temp = ''
        self.filename_save = pathdir + 'stats.csv'
        self.vol_stats_file = pathdir + 'emulated-historical-stats.csv'
        self.d_p_to_v_e_s_gap: dict[tuple[int,float],list[tuple[float,int,int,int]]] = {}
        self.d_3_gap_val = []
        self.d_3_geq_gap_hits = []
        self.d_3_geq_gap_shots = []
        self.d_3_shots = None
        self.d_5_gap_val = []
        self.d_5_geq_gap_hits = []
        self.d_5_geq_gap_shots = []
        self.d_5_shots = None
        pass

    def read_through_custom_counts(self):
        if path.isfile(self.filename_save) == False:
            raise ValueError('no result found')
        self.stats = sinter.read_stats_from_csv_files(self.filename_save)
        for stat in self.stats:
            if stat.decoder == 'desaturation':
                if stat.json_metadata['c'] == 'end2end-inplace-distillation' \
                    and stat.json_metadata['d1'] == 3 \
                    and stat.json_metadata['p'] == 0.001:
                    self.d_3_shots = stat.shots
                    custom = stat.custom_counts
                    gap_vals = set()
                    for key in custom:
                        gap_vals.add(int(key[1:]))
                    max_gap = max(gap_vals)
                    self.d_3_gap_val = [i for i in range(max_gap+1)]
                    for gap in self.d_3_gap_val:
                        try:
                            custom['C'+str(gap)]
                        except:
                            custom['C'+str(gap)] = 0
                        try:
                            custom['E'+str(gap)]
                        except:
                            custom['E'+str(gap)] = 0
                    e_counts = 0
                    c_counts = 0
                    for gap in reversed(self.d_3_gap_val):
                        e_counts += custom['E'+str(gap)]
                        c_counts += custom['C'+str(gap)]
                        self.d_3_geq_gap_hits.insert(0,e_counts)
                        self.d_3_geq_gap_shots.insert(0,e_counts+c_counts)
                elif stat.json_metadata['c'] == 'end2end-inplace-distillation' \
                    and stat.json_metadata['d1'] == 5 \
                    and stat.json_metadata['p'] == 0.001:
                    self.d_5_shots = stat.shots
                    custom = stat.custom_counts
                    gap_vals = set()
                    for key in custom:
                        gap_vals.add(int(key[1:]))
                    max_gap = max(gap_vals)
                    self.d_5_gap_val = [i for i in range(max_gap+1)]
                    for gap in self.d_5_gap_val:
                        try:
                            custom['C'+str(gap)]
                        except:
                            custom['C'+str(gap)] = 0
                        try:
                            custom['E'+str(gap)]
                        except:
                            custom['E'+str(gap)] = 0
                    e_counts = 0
                    c_counts = 0
                    for gap in reversed(self.d_5_gap_val):
                        e_counts += custom['E'+str(gap)]
                        c_counts += custom['C'+str(gap)]
                        self.d_5_geq_gap_hits.insert(0,e_counts)
                        self.d_5_geq_gap_shots.insert(0,e_counts+c_counts)
                    




    def read_through_vol_stats(self) -> None:
        self.vol_stats = sinter.read_stats_from_csv_files(self.vol_stats_file)
        for stat in self.vol_stats:
            if stat.decoder == 'desaturation':
                if stat.json_metadata['c'] == '2024 This Work (d1=3)':
                    try:
                        self.d_p_to_v_e_s_gap[(3,stat.json_metadata['p'])].append(
                            (stat.json_metadata['v'],
                             stat.errors,
                             stat.shots - stat.discards,
                             stat.json_metadata['gap'])
                        )
                    except:
                        self.d_p_to_v_e_s_gap[(3,stat.json_metadata['p'])] = [
                            (stat.json_metadata['v'],
                             stat.errors,
                             stat.shots - stat.discards,
                             stat.json_metadata['gap'])
                        ]
                elif stat.json_metadata['c'] == '2024 This Work (d1=5)':
                    try:
                        self.d_p_to_v_e_s_gap[(5,stat.json_metadata['p'])].append(
                            (stat.json_metadata['v'],
                             stat.errors,
                             stat.shots - stat.discards,
                             stat.json_metadata['gap'])
                        )
                    except:
                        self.d_p_to_v_e_s_gap[(5,stat.json_metadata['p'])] = [
                            (stat.json_metadata['v'],
                             stat.errors,
                             stat.shots - stat.discards,
                             stat.json_metadata['gap'])
                        ]

        pass


class SO_2d():

    def __init__(self,
                 file_name_no_append: str,
                 pathdir: str = './sample_results/'):
        self.filename_temp = pathdir + file_name_no_append +'_temp.csv'
        self.filename_save = pathdir + file_name_no_append + '_combined.csv'
        
        # if (path.isfile(self.filename_temp) or path.isfile(self.filename_save)) == False:
        #     raise ValueError('No results present.')


        self.shots = 0
        self.discards = 0
        self.errors = 0
        self.gaps_ec_2_counts: dict[tuple[int,int,bool],int] = {}
        self.gap_mono_vals = set()
        self.gap_vals = set()


        pass

    def clear_cache(self) -> None:
        if path.isfile(self.filename_temp):
            with open(self.filename_temp,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')


    def update(self) -> None:
        if path.isfile(self.filename_temp) == False:
            raise ValueError('no new data present')
        stats_temp = sinter.read_stats_from_csv_files(self.filename_temp)
        if len(stats_temp) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        stats_temp = stats_temp[0]
        if path.isfile(self.filename_save) == False:
            with open(self.filename_save,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')
                file.write(stats_temp.to_csv_line())
                file.write('\n')
        else:
            existing_stats = sinter.read_stats_from_csv_files(self.filename_save)
            if len(existing_stats) != 1:
                raise ValueError('More than one task result is stored. Unsupported behavior for now.')
            existing_stats = existing_stats[0]
            updated_stats = stats_temp + existing_stats
            with open(self.filename_save,'w') as file:
                file.write(sinter.CSV_HEADER)
                file.write('\n')
                file.write(updated_stats.to_csv_line())
                file.write('\n')
        
        with open(self.filename_temp,'w') as file:
            file.write(sinter.CSV_HEADER)
            file.write('\n')
    



    def read_through_custom_counts(self) -> None:
        if path.isfile(self.filename_save) == False:
            raise ValueError('no result found')
        self.stats = sinter.read_stats_from_csv_files(self.filename_save)
        if len(self.stats) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        self.stats = self.stats[0]
        self.shots = self.stats.shots
        self.discards = self.stats.discards
        self.errors = self.stats.errors
        self.c_at_mono_geq_gap = [] # 2d list
        self.e_at_mono_geq_gap = [] # 2d list
        # self.c_geq_mono_and_gap = []
        # self.e_geq_mono_and_gap = []
        if len(self.stats.custom_counts) > 0:
            for key in self.stats.custom_counts:
                if key[0] == 'C':
                    err = False
                elif key[0] == 'E':
                    err = True
                sep_index = key.index('|')
                gap_mono = int(key[1:sep_index])
                gap = int(key[sep_index+1:])
                counts = self.stats.custom_counts[key]
                self.gap_mono_vals.add(gap_mono)
                self.gap_vals.add(gap)
                self.gaps_ec_2_counts[(gap_mono,gap,err)] = counts
            
            gap_mono_min = 0
            gap_mono_max = max(self.gap_mono_vals)
            gap_min = 0 
            gap_max = max(self.gap_vals)
            self.gap_mono_list = [i for i in range(gap_mono_min,gap_mono_max+1)]
            self.gap_list = [i for i in range(gap_min,gap_max+1)]
            for gap_mono in self.gap_mono_list:
                for gap in self.gap_list:
                    try:
                        self.gaps_ec_2_counts[(gap_mono,gap,True)]
                    except:
                        self.gaps_ec_2_counts[(gap_mono,gap,True)] = 0
                    try:
                        self.gaps_ec_2_counts[(gap_mono,gap,False)]
                    except:
                        self.gaps_ec_2_counts[(gap_mono,gap,False)] = 0
            
            for gap_mono in self.gap_mono_list:
                c_row = []
                e_row = []
                c_counts = 0
                e_counts = 0
                for gap in reversed(self.gap_list):
                    c_counts += self.gaps_ec_2_counts[(gap_mono,gap,False)]
                    e_counts += self.gaps_ec_2_counts[(gap_mono,gap,True)]
                    c_row.insert(0,c_counts)
                    e_row.insert(0,e_counts)
                self.c_at_mono_geq_gap.append(c_row)
                self.e_at_mono_geq_gap.append(e_row)
            
            c_geq_gap_and_mono = []
            e_geq_gap_and_mono = []
            for gap in self.gap_list:
                c_row = []
                e_row = []
                c_counts = 0
                e_counts = 0
                for gap_mono in reversed(self.gap_mono_list):
                    c_counts += self.c_at_mono_geq_gap[gap_mono][gap]
                    e_counts += self.e_at_mono_geq_gap[gap_mono][gap]
                    c_row.insert(0,c_counts)
                    e_row.insert(0,e_counts)
                c_geq_gap_and_mono.append(c_row)
                e_geq_gap_and_mono.append(e_row)
            
            self.c_geq_mono_and_gap_arr = np.transpose(np.array(c_geq_gap_and_mono))
            self.e_geq_mono_and_gap_arr = np.transpose(np.array(e_geq_gap_and_mono))
                



    def rate(self, hits_list: list[int], shots_list: list[int]) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=1000).best)
        return ret_list

    def rate_ceil_binfit(self, hits_list: list[int],
                         shots_list: list[int],
                         param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).high)
        return ret_list

    def rate_floor_binfit(self, hits_list: list[int],
                          shots_list: list[int],
                          param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).low)
        return ret_list