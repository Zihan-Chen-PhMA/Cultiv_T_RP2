from csv_so_processor import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from volume_calc import *

gap_handler = Gap()
gap_handler.read_through_custom_counts()
gap_handler.read_through_vol_stats()


d_rp3 = 3
d_sc = 7
rp3_name = 'rp' + str(d_rp3)
sc_name = 'sc' + str(d_sc)
n_rounds_d3 = d_rp3 

filename_d3_e2e = 'rp_' + str(d_rp3) + '_sc_' + str(d_sc) + '_end2end_' \
                + str(n_rounds_d3+1) + '_full_rds' + '_so2d'
d3_end2end_res = SO_2d(filename_d3_e2e)
d3_end2end_res.read_through_custom_counts()

filename_d3_circ = 'rp_3_sc_7_end2end_4_full_rds'
circ_path = './circuit_garage/' + filename_d3_circ + '.stim'
discard_res_path = './sample_results/' + 'rp_3_T_cult' + '_combined.csv'


d3_vol_helper = Circ_rp3_T_end2end(circ_path,d_sc,n_rounds_d3)
d3_vol_helper.config_stages()
d3_vol_helper.active_qubits_calc()
d3_vol_helper.load_from_discard_tests(discard_res_path)


d_rp3 = 3
d_rp5 = 5
d_sc = 11
rp3_name = 'rp' + str(d_rp3)
rp5_name = 'rp' + str(d_rp5)
sc_name = 'sc' + str(d_sc)
n_rounds_d5 = d_rp5 
filename_d5_e2e = 'rp_' + str(d_rp3) \
                + '_rp_' + str(d_rp5) \
                + '_sc_' + str(d_sc) + '_end2end_' \
                + str(n_rounds_d5+1) + '_full_rds' + '_so2d'

d5_end2end_res = SO_2d(filename_d5_e2e)
d5_end2end_res.read_through_custom_counts()



filename_d5_circ = 'rp_3_rp_5_sc_11_end2end_6_full_rds'
circ_path_d5 = './circuit_garage/' + filename_d5_circ + '.stim'
discard_res_path_d5 = './sample_results/' + 'rp_3_rp_5_T_cult' + '_combined.csv'


d5_vol_helper = Circ_rp5_T_end2end(circ_path_d5,d_sc,n_rounds_d5)
d5_vol_helper.config_stages()
d5_vol_helper.active_qubits_calc()
d5_vol_helper.load_from_discard_tests(discard_res_path_d5)



d3_retry_arr_2d = 1/((d3_end2end_res.c_geq_mono_and_gap_arr 
                + d3_end2end_res.e_geq_mono_and_gap_arr)/d3_end2end_res.shots)
d3_err_rates_arr_2d = d3_end2end_res.e_geq_mono_and_gap_arr/(
                        d3_end2end_res.c_geq_mono_and_gap_arr + 
                        d3_end2end_res.e_geq_mono_and_gap_arr)
d3_hits_1d = np.reshape(d3_end2end_res.e_geq_mono_and_gap_arr,(-1)).astype(np.int64)
d3_shots_1d = np.reshape(d3_end2end_res.e_geq_mono_and_gap_arr 
                      + d3_end2end_res.c_geq_mono_and_gap_arr,(-1)).astype(np.int64)

# get best trade-off indices:
d3_retry_arr_100 = np.round((100*np.reshape(d3_retry_arr_2d,(-1)))).astype(np.int64)

d3_retry_arr_1d = np.reshape(d3_retry_arr_2d,(-1))
d3_err_rates_arr_1d = np.reshape(d3_err_rates_arr_2d,(-1))

d3_best_indices = []
d3_retry_2_err_rates: dict[int,list[tuple[int,float]]] = {}
for index in range(len(d3_retry_arr_100)):
    try:
        d3_retry_2_err_rates[d3_retry_arr_100[index]].append((index,d3_err_rates_arr_1d[index]))
    except:
        d3_retry_2_err_rates[d3_retry_arr_100[index]] = [(index,d3_err_rates_arr_1d[index])]

d3_retry_rates_100 = sorted(d3_retry_2_err_rates.keys())
for rate_100 in d3_retry_rates_100:
    index_temp = d3_retry_2_err_rates[rate_100][0][0]
    val_temp = d3_retry_2_err_rates[rate_100][0][1]
    for ind_err in d3_retry_2_err_rates[rate_100]:
        if ind_err[1] < val_temp:
            index_temp = ind_err[0]
            val_temp = ind_err[1]
    if d3_err_rates_arr_1d[index_temp] > 0 :
        d3_best_indices.append(index_temp)

d3_best_indices_distilled = [d3_best_indices[0]]
for index in d3_best_indices:
    if d3_err_rates_arr_1d[index] < d3_err_rates_arr_1d[d3_best_indices_distilled[-1]]:
        d3_best_indices_distilled.append(index)
d3_best_indices = []
d3_best_indices = [index for index in d3_best_indices_distilled]


d3_retry_rate_best = [d3_retry_arr_1d[index] for index in d3_best_indices]
d3_err_rate_best = [d3_err_rates_arr_1d[index] for index in d3_best_indices]

d3_hits_best = [d3_hits_1d[index] for index in d3_best_indices]
d3_shots_best = [d3_shots_1d[index] for index in d3_best_indices]
d3_err_rate_high = d3_end2end_res.rate_ceil_binfit(d3_hits_best,d3_shots_best)
d3_err_rate_low = d3_end2end_res.rate_floor_binfit(d3_hits_best,d3_shots_best)



d5_retry_arr_2d = 1/((d5_end2end_res.c_geq_mono_and_gap_arr 
                + d5_end2end_res.e_geq_mono_and_gap_arr)/d5_end2end_res.shots)
d5_err_rates_arr_2d = d5_end2end_res.e_geq_mono_and_gap_arr/(
                        d5_end2end_res.c_geq_mono_and_gap_arr + 
                        d5_end2end_res.e_geq_mono_and_gap_arr)
d5_hits_1d = np.reshape(d5_end2end_res.e_geq_mono_and_gap_arr,(-1)).astype(np.int64)
d5_shots_1d = np.reshape(d5_end2end_res.e_geq_mono_and_gap_arr 
                      + d5_end2end_res.c_geq_mono_and_gap_arr,(-1)).astype(np.int64)


# get best trade-off indices:
d5_retry_arr_100 = np.round((100*np.reshape(d5_retry_arr_2d,(-1)))).astype(np.int64)

d5_retry_arr_1d = np.reshape(d5_retry_arr_2d,(-1))
d5_err_rates_arr_1d = np.reshape(d5_err_rates_arr_2d,(-1))

d5_best_indices = []
d5_retry_2_err_rates: dict[int,list[tuple[int,float]]] = {}
for index in range(len(d5_retry_arr_100)):
    try:
        d5_retry_2_err_rates[d5_retry_arr_100[index]].append((index,d5_err_rates_arr_1d[index]))
    except:
        d5_retry_2_err_rates[d5_retry_arr_100[index]] = [(index,d5_err_rates_arr_1d[index])]

d5_retry_rates_100 = sorted(d5_retry_2_err_rates.keys())
for rate_100 in d5_retry_rates_100:
    index_temp = d5_retry_2_err_rates[rate_100][0][0]
    val_temp = d5_retry_2_err_rates[rate_100][0][1]
    for ind_err in d5_retry_2_err_rates[rate_100]:
        if ind_err[1] < val_temp:
            index_temp = ind_err[0]
            val_temp = ind_err[1]
    if d5_err_rates_arr_1d[index_temp] > 0:
        d5_best_indices.append(index_temp)

# further distillation:
d5_best_indices_distilled = [d5_best_indices[0]]
for index in d5_best_indices:
    if d5_err_rates_arr_1d[index] < d5_err_rates_arr_1d[d5_best_indices_distilled[-1]]:
        d5_best_indices_distilled.append(index)

d5_best_indices = []
d5_best_indices = [index for index in d5_best_indices_distilled]


d5_retry_rate_best = [d5_retry_arr_1d[index] for index in d5_best_indices]
d5_err_rate_best = [d5_err_rates_arr_1d[index] for index in d5_best_indices]
d5_hits_best = [d5_hits_1d[index] for index in d5_best_indices]
d5_shots_best = [d5_shots_1d[index] for index in d5_best_indices]
d5_err_rate_high = d5_end2end_res.rate_ceil_binfit(d5_hits_best,d5_shots_best)
d5_err_rate_low = d5_end2end_res.rate_floor_binfit(d5_hits_best,d5_shots_best)


msc_color = '#053691'
rp_color = '#ff8070'
# rp_color = '#f55e19'






d3_selected_indices = []
for i, err_rate in zip(range(len(d3_err_rate_best)),d3_err_rate_best):
    if len(d3_selected_indices) == 0:
        d3_selected_indices.append(i)
    else:
        if err_rate < 2e-6:
            if err_rate/d3_err_rate_best[d3_selected_indices[-1]] > 0.9:
                pass
            else:
                d3_selected_indices.append(i)
            continue
        if err_rate/d3_err_rate_best[d3_selected_indices[-1]] > 0.9:
            pass
        else:
            d3_selected_indices.append(i)

selected_d3_retry_rate = [d3_retry_rate_best[i] for i in d3_selected_indices]
selected_d3_err_rate_best = [d3_err_rate_best[i] for i in d3_selected_indices]
selected_d3_err_rate_high = [d3_err_rate_high[i] for i in d3_selected_indices]
selected_d3_err_rate_low = [d3_err_rate_low[i] for i in d3_selected_indices]


d5_selected_indices = []
for i, err_rate in zip(range(len(d5_err_rate_best)),d5_err_rate_best):
    if err_rate < 2.5e-10:
        continue
    else:
        d5_selected_indices.append(i)
        

selected_d5_retry_rate = [d5_retry_rate_best[i] for i in d5_selected_indices]
selected_d5_err_rate_best = [d5_err_rate_best[i] for i in d5_selected_indices]
selected_d5_err_rate_high = [d5_err_rate_high[i] for i in d5_selected_indices]
selected_d5_err_rate_low = [d5_err_rate_low[i] for i in d5_selected_indices]



size = (6.4,4.8)


fig, ax = plt.subplots()
fig.set_size_inches(size)

ax.set_axisbelow(True)
ax.yaxis.grid(True,which='major',color='gray')
ax.yaxis.grid(True,which='minor',color='#cccccc')
ax.xaxis.grid(True,which='major',color='gray')
ax.plot(selected_d3_retry_rate,selected_d3_err_rate_best,marker='^',
        linestyle='dotted',
        color=rp_color)

ax.fill_between(selected_d3_retry_rate,selected_d3_err_rate_low,selected_d3_err_rate_high,
                 color=rp_color,alpha=0.6)

msc_d3_gap_retry_rates = []
msc_d3_gap_rate_best = []
msc_d3_gap_rate_high = []
msc_d3_gap_rate_low = []
msc_d5_gap_rate_best = []
msc_d5_gap_rate_high = []
msc_d5_gap_rate_low = []
msc_d5_gap_retry_rates = []

for hits, shots in zip(gap_handler.d_3_geq_gap_hits,
                       gap_handler.d_3_geq_gap_shots):
    msc_d3_gap_retry_rates.append(gap_handler.d_3_shots/shots)
    rate_fit = sinter.fit_binomial(num_hits=hits,num_shots=shots,max_likelihood_factor=1000)
    msc_d3_gap_rate_best.append(rate_fit.best)
    msc_d3_gap_rate_high.append(rate_fit.high)
    msc_d3_gap_rate_low.append(rate_fit.low)

for hits, shots in zip(gap_handler.d_5_geq_gap_hits,
                       gap_handler.d_5_geq_gap_shots):
    msc_d5_gap_retry_rates.append(gap_handler.d_5_shots/shots)
    rate_fit = sinter.fit_binomial(num_hits=hits,num_shots=shots,max_likelihood_factor=1000)
    msc_d5_gap_rate_best.append(rate_fit.best)
    msc_d5_gap_rate_high.append(rate_fit.high)
    msc_d5_gap_rate_low.append(rate_fit.low)


# ax.scatter(msc_d3_gap_retry_rates,msc_d3_gap_rate_best,marker='h')
ax.plot(msc_d3_gap_retry_rates, msc_d3_gap_rate_best, marker='h',
        linestyle = 'dotted',
        color=msc_color)
ax.fill_between(msc_d3_gap_retry_rates, msc_d3_gap_rate_low,
                msc_d3_gap_rate_high, color=msc_color, alpha=0.4)


ax.plot(selected_d5_retry_rate, selected_d5_err_rate_best, marker='v',
        linestyle='dotted', color=rp_color)
ax.fill_between(selected_d5_retry_rate, selected_d5_err_rate_low,
                selected_d5_err_rate_high, color=rp_color, alpha=0.6)

# ax.scatter(d5_retry_rate_best,d5_err_rate_best)
# ax.fill_between(d5_retry_rate_best,d5_err_rate_low,d5_err_rate_high,
#                  color='orange',alpha=0.4)


ax.plot(msc_d5_gap_retry_rates,msc_d5_gap_rate_best,marker='H',
        linestyle='dotted',color=msc_color)
ax.fill_between(msc_d5_gap_retry_rates,msc_d5_gap_rate_low,msc_d5_gap_rate_high,
                color = msc_color, alpha=0.4)

ax.set_xlabel('Expected Attempts per Kept Shot',fontsize=12)
ax.set_ylabel('Logilca Error Rate',fontsize=12)
ax.tick_params(labelsize=12)
ax.set_yscale('logit')
ax.set_xscale('log')
fig.savefig('LER vs expected attempts.pdf',bbox_inches='tight')




for stage in d5_vol_helper.stages:
    print(stage.stage_name, stage.active_qubits, stage.survival_rate, stage.det_stage)

# print(vol_helper.volume_calc(1/2.3),vol_helper.volume_calc(1/2.5))





for stage in d3_vol_helper.stages:
    print(stage.stage_name, stage.active_qubits, stage.survival_rate, stage.det_stage)










msc_d3_v = []    
msc_d3_errors_v = []
msc_d3_shots_v = []
msc_d3_rate_v = []
msc_d3_rate_high_v = []
msc_d3_rate_low_v = []
for item in gap_handler.d_p_to_v_e_s_gap[(3,0.001)]:
    msc_d3_v.append(item[0])
    msc_d3_errors_v.append(round(item[1]/2))
    msc_d3_shots_v.append(item[2])
    err_rate_fit = sinter.fit_binomial(num_hits=round(item[1]/2),
                                       num_shots=item[2],max_likelihood_factor=1000)
    msc_d3_rate_v.append(err_rate_fit.best)
    msc_d3_rate_high_v.append(err_rate_fit.high)
    msc_d3_rate_low_v.append(err_rate_fit.low)

msc_d5_v = []    
msc_d5_errors_v = []
msc_d5_shots_v = []
msc_d5_rate_v = []
msc_d5_rate_high_v = []
msc_d5_rate_low_v = []
for item in gap_handler.d_p_to_v_e_s_gap[(5,0.001)]:
    msc_d5_v.append(item[0])
    msc_d5_errors_v.append(round(item[1]/2))
    msc_d5_shots_v.append(item[2])
    err_rate_fit = sinter.fit_binomial(num_hits=round(item[1]/2),
                                       num_shots=item[2],max_likelihood_factor=1000)
    msc_d5_rate_v.append(err_rate_fit.best)
    msc_d5_rate_high_v.append(err_rate_fit.high)
    msc_d5_rate_low_v.append(err_rate_fit.low)



d3_selected_indices = []
for i, err_rate in zip(range(len(d3_err_rate_best)),d3_err_rate_best):
    if len(d3_selected_indices) == 0:
        d3_selected_indices.append(i)
    else:
        if err_rate < 2e-6:
            if err_rate/d3_err_rate_best[d3_selected_indices[-1]] > 0.87:
                pass
            else:
                d3_selected_indices.append(i)
            continue
        if err_rate/d3_err_rate_best[d3_selected_indices[-1]] > 0.75:
            pass
        else:
            d3_selected_indices.append(i)

selected_d3_retry_rate = [d3_retry_rate_best[i] for i in d3_selected_indices]
selected_d3_err_rate_best = [d3_err_rate_best[i] for i in d3_selected_indices]
selected_d3_err_rate_high = [d3_err_rate_high[i] for i in d3_selected_indices]
selected_d3_err_rate_low = [d3_err_rate_low[i] for i in d3_selected_indices]




d5_selected_indices = []
for i, err_rate in zip(range(len(d5_err_rate_best)),d5_err_rate_best):
    if err_rate > 2e-6:
        continue
    elif len(d5_selected_indices) == 0:
        d5_selected_indices.append(i)
    else:
        if 9e-10 < err_rate < 1e-9:
            d5_selected_indices.append(i)
            continue
        if err_rate < 3e-10:
            continue
        if err_rate < 5e-10:
            if err_rate/d5_err_rate_best[d5_selected_indices[-1]] < 0.9:
                d5_selected_indices.append(i)
                continue
        if err_rate < 4e-9:
            if err_rate/d5_err_rate_best[d5_selected_indices[-1]] > 0.8:
                pass
            elif err_rate < 5e-10:
                d5_selected_indices.append(i)
            else:
                d5_selected_indices.append(i)
            continue
        if err_rate/d5_err_rate_best[d5_selected_indices[-1]] > 0.7:
            pass
        else:
            d5_selected_indices.append(i)

selected_d5_retry_rate = [d5_retry_rate_best[i] for i in d5_selected_indices]
selected_d5_err_rate_best = [d5_err_rate_best[i] for i in d5_selected_indices]
selected_d5_err_rate_high = [d5_err_rate_high[i] for i in d5_selected_indices]
selected_d5_err_rate_low = [d5_err_rate_low[i] for i in d5_selected_indices]




fig, ax = plt.subplots()
fig.set_size_inches(size)

ax.set_axisbelow(True)
ax.yaxis.grid(True,which='major',color='gray')
ax.yaxis.grid(True,which='minor',color='#cccccc')
ax.xaxis.grid(True,which='major',color='gray')
d3_vol = [d3_vol_helper.volume_calc(1/retry_rate) for retry_rate in selected_d3_retry_rate]
d3_errbar = [[],[]]
for err, err_high, err_low in zip(selected_d3_err_rate_best,
                                  selected_d3_err_rate_high,
                                  selected_d3_err_rate_low):
    d3_errbar[0].append(err-err_low)
    d3_errbar[1].append(err_high-err)

ax.errorbar(d3_vol,selected_d3_err_rate_best,yerr=d3_errbar,ls='None',color=rp_color,
            alpha=0.6)

# ax.fill_between(d5_vol,selected_d5_err_rate_low,selected_d5_err_rate_high,
#                  color=rp_color,alpha=0.4)
ax.scatter(d3_vol,selected_d3_err_rate_best,marker='^',color=rp_color,s=30)

# ax.fill_between(d3_vol,d3_err_rate_low,d3_err_rate_high,
#                  color=rp_color,alpha=0.4)
# ax.scatter(d3_vol,d3_err_rate_best,marker='o',s=20,color=rp_color)


msc_rate_errbar_d3 = [[],[]]
for err, err_high, err_low in zip(msc_d3_rate_v, 
                                  msc_d3_rate_high_v,
                                  msc_d3_rate_low_v):
    msc_rate_errbar_d3[0].append(err-err_low)
    msc_rate_errbar_d3[1].append(err_high-err)

ax.errorbar(msc_d3_v,y=msc_d3_rate_v, yerr=msc_rate_errbar_d3,
            fmt='h',color=msc_color)

d5_vol = [d5_vol_helper.volume_calc(1/retry_rate) for retry_rate in selected_d5_retry_rate]
d5_errbar = [[],[]]
for err, err_high, err_low in zip(selected_d5_err_rate_best,
                                  selected_d5_err_rate_high,
                                  selected_d5_err_rate_low):
    d5_errbar[0].append(err-err_low)
    d5_errbar[1].append(err_high-err)
ax.errorbar(d5_vol,selected_d5_err_rate_best,yerr=d5_errbar,ls='None',color=rp_color,alpha=0.6)

# ax.fill_between(d5_vol,selected_d5_err_rate_low,selected_d5_err_rate_high,
#                  color=rp_color,alpha=0.4)
ax.scatter(d5_vol,selected_d5_err_rate_best,marker='v',color=rp_color,s=30)



msc_rate_errbar_d5 = [[],[]]
for err, err_high, err_low in zip(msc_d5_rate_v, 
                                  msc_d5_rate_high_v,
                                  msc_d5_rate_low_v):
    msc_rate_errbar_d5[0].append(err-err_low)
    msc_rate_errbar_d5[1].append(err_high-err)

ax.errorbar(msc_d5_v,y=msc_d5_rate_v, yerr=msc_rate_errbar_d5,
            linestyle='None',color=msc_color,alpha=0.4)
ax.scatter(msc_d5_v,msc_d5_rate_v,marker='H',color=msc_color)




ax.set_xlabel('Expected Space-time Volume',fontsize=12)
ax.set_ylabel('Logical Error Rate',fontsize=12)
ax.tick_params(labelsize=12)
ax.set_yscale('logit')
ax.set_xscale('log')

fig.savefig('LER vs expected spacetime volume.pdf',bbox_inches='tight')


def fmt(x):
    return '{:.1e}'.format(x)


fig, ax = plt.subplots()
fig.set_size_inches(size)

d3_gap_mono_max = max(d3_end2end_res.gap_mono_vals)
d3_gap_max = max(d3_end2end_res.gap_vals)

d3_so_indices = [d3_best_indices[i] for i in d3_selected_indices]
selected_d3_gap_mono = [index // (d3_gap_max+1) for index in d3_so_indices]
selected_d3_gap = [index % (d3_gap_max+1) for index in d3_so_indices]

y = np.arange(d3_gap_mono_max+1)
x = np.arange(d3_gap_max+1)

x_grid, y_grid = np.meshgrid(x,y)
norm = colors.LogNorm(np.min(d3_err_rates_arr_2d),np.max(d3_err_rates_arr_2d))
ax.pcolormesh(x_grid,y_grid, d3_err_rates_arr_2d,
              norm=norm, cmap='RdYlGn_r',shading='nearest')

levels = np.array([1.5e-6])
CS = ax.contour(x_grid,y_grid,d3_err_rates_arr_2d, levels=levels, colors='k')
ax.clabel(CS,CS.levels,fontsize=12,fmt=fmt)

ax.scatter(selected_d3_gap,selected_d3_gap_mono,facecolors='none',edgecolors=rp_color)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdYlGn_r'),ax=ax)
cbar.set_label('Logical Error Rate',fontsize=12)
cbar.ax.tick_params(labelsize=12)
ax.set_xlabel(r'$\phi_{bd}$ (dB)',fontsize=12)
ax.set_ylabel(r'$\phi_{\mathbb{RP}^2}$ (dB)',fontsize=12)
ax.tick_params(labelsize=12)
fig.savefig('soft_outputs_vs_LER_d3.pdf',bbox_inches='tight')


fig, ax = plt.subplots()
fig.set_size_inches(size)

d3_gap_mono_max = max(d3_end2end_res.gap_mono_vals)
d3_gap_max = max(d3_end2end_res.gap_vals)
y = np.arange(d3_gap_mono_max+1)
x = np.arange(d3_gap_max+1)

x_grid, y_grid = np.meshgrid(x,y)
norm = colors.LogNorm(np.min(d3_retry_arr_2d),np.max(d3_retry_arr_2d))
ax.pcolormesh(x_grid,y_grid, d3_retry_arr_2d,
              norm=norm, cmap='Spectral_r',shading='nearest')
levels = np.array([2.4])
CS = ax.contour(x_grid,y_grid,d3_retry_arr_2d, levels=levels, colors='k')
ax.clabel(CS,CS.levels,fontsize=12)

ax.scatter(selected_d3_gap,selected_d3_gap_mono,facecolors='none',edgecolors=rp_color)

ax.set_xlabel(r'$\phi_{bd}$ (dB)',fontsize=12)
ax.set_ylabel(r'$\phi_{\mathbb{RP}^2}$ (dB)',fontsize=12)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap='Spectral_r'),ax=ax)
cbar.set_label('Expected Attempts per Kept Shot',fontsize=12)
cbar.ax.tick_params(which='both',labelsize=12)
ax.tick_params(labelsize=12)
fig.savefig('soft_outputs_vs_spacetime_cost_d3.pdf',bbox_inches='tight')




fig, ax = plt.subplots()
fig.set_size_inches(size)

d5_gap_mono_max = max(d5_end2end_res.gap_mono_vals)
d5_gap_max = max(d5_end2end_res.gap_vals)

d5_so_indices = [d5_best_indices[i] for i in d5_selected_indices]
selected_d5_gap_mono = [index // (d5_gap_max+1) for index in d5_so_indices]
selected_d5_gap = [index % (d5_gap_max+1) for index in d5_so_indices]

y = np.arange(d5_gap_mono_max+1)
x = np.arange(d5_gap_max+1)

x_grid, y_grid = np.meshgrid(x,y)
norm = colors.LogNorm(np.min(d5_err_rates_arr_2d),np.max(d5_err_rates_arr_2d))
ax.pcolormesh(x_grid,y_grid, d5_err_rates_arr_2d,
              norm=norm, cmap='RdYlGn_r',shading='nearest')
levels = np.array([1e-9])
CS = ax.contour(x_grid,y_grid,d5_err_rates_arr_2d, levels=levels, colors='k')
ax.clabel(CS,CS.levels,fontsize=12,fmt=fmt)

ax.scatter(selected_d5_gap,selected_d5_gap_mono,facecolors='none',edgecolors=rp_color)

cbar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdYlGn_r'),ax=ax)
cbar.set_label('Logical Error Rate',fontsize=12)
cbar.ax.tick_params(labelsize=12)

ax.set_xlabel(r'$\phi_{bd}$ (dB)',fontsize=12)
ax.set_ylabel(r'$\phi_{\mathbb{RP}^2}$ (dB)',fontsize=12)
ax.tick_params(labelsize=12)
fig.savefig('soft_outputs_vs_LER_d5.pdf',bbox_inches='tight')



fig, ax = plt.subplots()
fig.set_size_inches(size)

d5_gap_mono_max = max(d5_end2end_res.gap_mono_vals)
d5_gap_max = max(d5_end2end_res.gap_vals)
y = np.arange(d5_gap_mono_max+1)
x = np.arange(d5_gap_max+1)

x_grid, y_grid = np.meshgrid(x,y)
norm = colors.LogNorm(np.min(d5_retry_arr_2d),np.max(d5_retry_arr_2d))
ax.pcolormesh(x_grid,y_grid, d5_retry_arr_2d,
              norm=norm, cmap='Spectral_r',shading='nearest')

levels = np.array([15],dtype=np.int64)
CS = ax.contour(x_grid,y_grid,d5_retry_arr_2d, levels=levels, colors='k')
ax.clabel(CS,CS.levels,fontsize=12)
ax.scatter(selected_d5_gap,selected_d5_gap_mono,facecolors='none',edgecolors=rp_color)

ax.set_xlabel(r'$\phi_{bd}$ (dB)',fontsize=12)
ax.set_ylabel(r'$\phi_{\mathbb{RP}^2}$ (dB)',fontsize=12)
ax.tick_params(labelsize=12)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap='Spectral_r'),ax=ax)
cbar.set_label('Expected Attempts per Kept Shot',fontsize=12)
cbar.ax.tick_params(which='both',labelsize=12)
fig.savefig('soft_outputs_vs_spacetime_cost_d5.pdf',bbox_inches='tight')


fig, ax = plt.subplots()
fig.set_size_inches(size)
qubits_color = '#e06700'
survival_color = '#00aa05'
alpha_qubit = 0.8
alpha=0.5

bins = []
widths = []
active_qubits = []
survival_rate = []
stage_names = []

bin_temp = 0

for stage in d3_vol_helper.stages:
    bins.append(bin_temp)
    if 'MorphBack&Grow' in stage.stage_name:
        bin_temp += 3
        widths.append(3)
    elif 'MorphTo' in stage.stage_name:
        bin_temp += 1
        widths.append(1)
    else:
        bin_temp += 2
        widths.append(2)
    
    active_qubits.append(stage.active_qubits)
    survival_rate.append(stage.survival_rate)

active_qubits[-1] = active_qubits[-2]


thres = 1.5e-6
ind_temp = 0
for i in range(len(selected_d3_err_rate_best)):
    if selected_d3_err_rate_best[i] < thres:
        print(selected_d3_err_rate_best[i])
        survival_rate[-1] = 1/selected_d3_retry_rate[i]
        break

for stage in d3_vol_helper.stages:
    if 'MorphBack&Grow' in stage.stage_name:
        if stage.active_qubits == max(active_qubits):
            stage_names.append('Morph &\nExpand')
        else:
            stage_names.append('Morph &\nGrow')
    elif 'MorphTo' in stage.stage_name:
        stage_names.append('Morph')
    elif 'Injection' in stage.stage_name:
        stage_names.append('Inject T')
    elif 'T_check' in stage.stage_name:
        stage_names.append('Check T')
    elif 'PerfMeas' in stage.stage_name:
        stage_names.append('Perfect\nMeas.')
    elif 'SE' in stage.stage_name:
        stage_names.append('SE')

# bins.append(bin_temp)

ax.set_ylabel('Activated Qubits',fontsize=12,color=qubits_color)
ax.tick_params(axis='y',color=qubits_color,
               labelcolor=qubits_color)
ax.bar(bins,active_qubits,widths,align='edge',
       facecolor=colors.to_rgba(qubits_color,alpha_qubit),
       edgecolor=qubits_color)
ax.set_xticks(bins + [bin_temp])
ax.tick_params(axis='x',labelbottom=False)
ax.set_xlim(left=0,right=bin_temp)
ax.set_ylim(bottom=0,top=max(active_qubits))

bin_centers = []
for i in range(len(bins)):
    if i+1 != len(bins):
        bin_centers.append((bins[i]+bins[i+1])/2)
    else:
        bin_centers.append((bins[-1]+bin_temp)/2)

for bin_center, name in zip(bin_centers,stage_names):
    ax.annotate(name,xy=(bin_center,0.05*max(active_qubits)),
                xytext=(bin_center,-0.01*max(active_qubits)),
                rotation=90, horizontalalignment='center',
                verticalalignment='top',fontsize=10)

ax2 = ax.twinx()
ax2.set_ylabel('Survival Rate',fontsize=12,color=survival_color)
ax2.tick_params(axis='y',color=survival_color,
                labelcolor=survival_color)
ax2.bar(bins,survival_rate,widths,align='edge',
        facecolor=colors.to_rgba(survival_color,alpha),
        edgecolor=survival_color)
ax2.set_ylim(bottom=0,top=1)
fig.savefig('MSC_3_process_analysis.pdf',bbox_inches='tight')





fig, ax = plt.subplots()
fig.set_size_inches(size)
qubits_color = '#e06700'
survival_color = '#00aa05'
alpha_qubit = 0.8
alpha=0.5

bins = []
widths = []
active_qubits = []
survival_rate = []
stage_names = []

bin_temp = 0

for stage in d5_vol_helper.stages:
    bins.append(bin_temp)
    if 'MorphBack&Grow' in stage.stage_name:
        bin_temp += 3
        widths.append(3)
    elif 'MorphTo' in stage.stage_name:
        bin_temp += 1
        widths.append(1)
    else:
        bin_temp += 2
        widths.append(2)
    
    active_qubits.append(stage.active_qubits)
    survival_rate.append(stage.survival_rate)

active_qubits[-1] = active_qubits[-2]

thres = 1e-9
ind_temp = 0
for i in range(len(selected_d5_err_rate_best)):
    if selected_d5_err_rate_best[i] < thres:
        print(selected_d5_err_rate_best[i])
        survival_rate[-1] = 1/selected_d5_retry_rate[i]
        break


for stage in d5_vol_helper.stages:
    if 'MorphBack&Grow' in stage.stage_name:
        if stage.active_qubits == max(active_qubits):
            stage_names.append('Morph &\nExpand')
        else:
            stage_names.append('Morph &\nGrow')
    elif 'MorphTo' in stage.stage_name:
        stage_names.append('Morph')
    elif 'Injection' in stage.stage_name:
        stage_names.append('Inject T')
    elif 'T_check' in stage.stage_name:
        stage_names.append('Check T')
    elif 'PerfMeas' in stage.stage_name:
        stage_names.append('Perfect\nMeas.')
    elif 'SE' in stage.stage_name:
        stage_names.append('SE')

# bins.append(bin_temp)

ax.set_ylabel('Activated Qubits',fontsize=12,color=qubits_color)
ax.tick_params(axis='y',color=qubits_color,
               labelcolor=qubits_color)
ax.bar(bins,active_qubits,widths,align='edge',
       facecolor=colors.to_rgba(qubits_color,alpha_qubit),
       edgecolor=qubits_color)
ax.set_xticks(bins + [bin_temp])
ax.tick_params(axis='x',labelbottom=False)
ax.set_xlim(left=0,right=bin_temp)
ax.set_ylim(bottom=0,top=max(active_qubits))

bin_centers = []
for i in range(len(bins)):
    if i+1 != len(bins):
        bin_centers.append((bins[i]+bins[i+1])/2)
    else:
        bin_centers.append((bins[-1]+bin_temp)/2)

for bin_center, name in zip(bin_centers,stage_names):
    ax.annotate(name,xy=(bin_center,0.05*max(active_qubits)),
                xytext=(bin_center,-0.01*max(active_qubits)),
                rotation=90, horizontalalignment='center',
                verticalalignment='top',fontsize=10)

ax2 = ax.twinx()
ax2.set_ylabel('Survival Rate',fontsize=12,color=survival_color)
ax2.tick_params(axis='y',color=survival_color,
                labelcolor=survival_color)
ax2.bar(bins,survival_rate,widths,align='edge',
        facecolor=colors.to_rgba(survival_color,alpha),
        edgecolor=survival_color)
ax2.set_ylim(bottom=0,top=1)
fig.savefig('MSC_5_process_analysis.pdf',bbox_inches='tight')
