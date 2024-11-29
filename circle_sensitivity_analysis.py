# -*- coding: utf-8 -*-
"""
Code for  "Spatiotemporal determinants of stretch-activated channel-induced 
re-entry in ventricular tissue: an in-silico study" 
by Melania Buonocunto, Aurore Lyon, Tammo Delhaas, Joost Lumens, and Jordi Heijman
contact: work: m.buonocunto@maastrichtuniversity.nl ; melania.buonocunto@medunigraz.at
personal: melania.buonoc@gmail.com  
"""

import numpy as np
import myokit
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pickle as pkl
import datetime
import time
import multiprocessing



def run_function_in_parallel(target_function, data_to_run, max_processes = 10): 
    """
    Run function in parallel

    Parameters
    ----------
    target_function : callable
        Function with arguments (i_process, data, q_inbetween, q_final)
    data_to_run : list of any
        Each item is send to the target function
    max_processes : int, optional
        Maximum number of processes. The default is 10.

    Returns
    -------
    results : list of any
        List of results generated by target_function

    """
    # create a multiprocessing Queue
    q_inbetween = multiprocessing.Queue()
    q_result = multiprocessing.Queue()

    # store processes in list
    processes = []

    # initiate iteration
    n_done = 0
    i_process = 0

    while True:
        # start a new process when <max_processes processes are active
        while len(processes) < max_processes and i_process < len(data_to_run):
            data = data_to_run[i_process]

            process = multiprocessing.Process(target=target_function, args=(i_process,data,q_inbetween, q_result))
            processes.append(process)
            process.start()

            i_process += 1

        # handle when process is finished
        for process in processes:
            if not process.is_alive():
                processes.remove(process)
                n_done +=1

        # Print if data in queue, otherwise sleep
        if not q_inbetween.empty():
            print(q_inbetween.get())
        elif n_done == len(data_to_run):
            break
        else:
            time.sleep(1)

    # Finished! Now collect data
    results = [None for _ in range(len(data_to_run))]
    for i in range(len(data_to_run)):
        data = q_result.get()
        results[data[0]] = data[1]
    return results

############################
#Multiparameter analysis

#total ncells
n_list = [500]
#ns/n ratio stretched cells/total ncells
ns_ratio_list = np.arange(280,510,10) 
#stretch time of application
ts_list =  np.arange(350,390,10) 


n = 0
ns_ratio = 0
ts = 0
tot = 0
df = pd.DataFrame(columns = ['n', 'ns_ratio', 'ts'])

for n in n_list:
    for ns_ratio in ns_ratio_list:
        for ts in ts_list:
            tot = tot+1 
            df_temp = pd.DataFrame({'n': n, 'ns_ratio': ns_ratio, 'ts': ts, 'count': tot}, index=[1])
            df = pd.concat([df_temp,df])

df = pd.concat([df])
df.to_csv('simulation_stretch_centre.csv', index=False)


if __name__ == '__main__':

    from circle_TORORdSACs_2Dsens import model
    results = run_function_in_parallel(model, df.values)
    df = pd.DataFrame(results)
    now = datetime.datetime.now()
    print(str(now))
    print(df)
    dfn = df.to_csv('simulation_stretch_centre.csv',index=False)
          