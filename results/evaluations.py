
# basics
import os
import re
import csv
import json
import matplotlib.pyplot as plt

# sklearn
from sklearn.metrics import auc


"""
    Utilities
"""
def isfloat(value):
    value = value.lstrip('-')
    return re.match("^\d+?\.\d+?$", value)

def read_csv_results(filename):
    ret_data = {
        'parameters': {},   # hyper-parameters
        'controls'  : {},   # control data: as { ess: [...], non-ess: [...] } format
        'results'   : {},   # result data : as { rxn: { ess: [...], non-ess: [...] }, ... }
    }
    # read
    with open(filename, 'rb') as inputfile:
        csv_reader = csv.reader(inputfile)
        rxn_remove = None
        for idx, each_line in enumerate(csv_reader):
            # parameters read
            if 0 <= idx < 6:
                param_value = each_line[1].strip()
                if isfloat(param_value): param_value = float(param_value)
                ret_data['parameters'][each_line[0]] = param_value
            # control data read
            if idx == 6: 
                data = each_line[1:]
                data = [each_data for each_data in data if each_data]
                ret_data['controls']['ess']     = data
            if idx == 7:
                data = each_line[1:]
                data = [each_data for each_data in data if each_data]
                ret_data['controls']['non-ess'] = data
            # results data read (from 9th row)
            if idx > 8:
                if idx % 3 == 0:
                    rxn_remove = each_line[0].strip()
                    rxn_remove = rxn_remove.replace('[', '')
                    rxn_remove = rxn_remove.replace(']', '')
                    ret_data['results'][rxn_remove] = {}
                if idx % 3 == 1:
                    data = each_line[1:]
                    data = [each_data for each_data in data if each_data]
                    ret_data['results'][rxn_remove]['ess']     = data
                if idx % 3 == 2:
                    data = each_line[1:]
                    data = [each_data for each_data in data if each_data]
                    ret_data['results'][rxn_remove]['non-ess'] = data
            # end read
    return ret_data

def compute_fprtprs(result):
    ret_data = {
        'parameters' : result['parameters'],
        'results'    : {},
    }
    # control sets
    control_ess  = set(result['controls']['ess'])
    control_ness = set(result['controls']['non-ess'])
    # compute the fpr and tprs
    for each_rxn, each_data in result['results'].iteritems():
        # compute TP, TN, FP, FN
        cur_ess  = set(each_data['ess'])
        cur_ness = set(each_data['non-ess'])
        true_positive  = float(len(control_ess.intersection(cur_ess)))  / len(control_ess)
        true_negative  = float(len(control_ness.intersection(cur_ness)))/ len(control_ness)
        false_positive = float(len(control_ness.intersection(cur_ess))) / len(control_ness)
        false_negative = float(len(control_ess.intersection(cur_ness))) / len(control_ess)
        # compute tpr and fpr
        cur_fpr = 0
        if (true_positive + false_negative): cur_fpr = true_positive  / (true_positive + false_negative)
        cur_tpr = 0
        if (false_positive + true_negative): cur_tpr = false_positive / (false_positive + true_negative)
        # store
        ret_data['results'][each_rxn] = [cur_fpr, cur_tpr]
    return ret_data

def extract_total_rxns(results):
    rxn_set = []
    for each_res in results:
        cur_res  = each_res['results']
        rxn_set += cur_res.keys()
    return set(rxn_set)

def aggregate_fprtprs_per_rxn(rxns, fprtprs):
    ret_data = { rxn: [] for rxn in rxns }
    # store the fpr and tprs
    for fprtpr in fprtprs:
        cur_res = fprtpr['results']
        for rxn, ft_rate in cur_res.iteritems():
            ret_data[rxn].append((ft_rate[0], ft_rate[1]))
    return ret_data

def draw_roc_curve(rxn, fprs, tprs, auc):
    # draw and save
    plt.plot(cur_fprs, cur_tprs, color='darkorange')
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('FP rates')
    plt.ylabel('TP rates')
    plt.title('The ROC curve for %s (AUC: %s)' % (rxn, cur_auc))
    plt.savefig('../evaluations/%s.eps' % (rxn))
    plt.clf()
    #Fin.

def save_to_file(filename, rxn_fptps):
    with open(filename, 'wb') as outputfile:
        for rxn, fptp in rxn_fptps.iteritems():
            output_str = ['(%03.6f, %03.6f)' % (data[0], data[1]) for data in fptp]
            outputfile.write('[%s] : \n' % (rxn))
            for each_str in output_str:
                outputfile.write('  %s\n' % (each_str))
    #Fin.


"""
    Evaluation codes
"""
def evaluations():
    print ' .. do evaluations of the model accuracy '

    # read the data files
    datafiles = [filename for filename in os.listdir('.') if filename.endswith('.csv')]
    print '    [%4s] total results files exist ' % (len(datafiles))

    # read the datafiles one by one
    results   = [read_csv_results(filename) for filename in datafiles]
    print '    [%4s] total results read from the csv files ' % (len(results))

    # compute the TP rate and FP rate for each data
    res_fptps = [compute_fprtprs(each_res) for each_res in results]
    print '    [%4s] total results about fpr/tpr for each rxn ' % (len(res_fptps))

    # compute the total reactions from results
    rxn_total = extract_total_rxns(results)
    print '    [%4s] total number of reactions ' % (len(rxn_total))

    # aggregate the fpr/tprs per rxn
    rxn_fptps = aggregate_fprtprs_per_rxn(rxn_total, res_fptps)
    print '    [%4s] total number of reactions and fpr/tprs ' % (len(rxn_fptps))

    # draw the ROC curve with computed the AUC in each rxn
    rxn_fprtpr= {}
    for rxn, fprtprs in rxn_fptps.iteritems():
        cur_fprs = [data[0] for data in fprtprs]
        cur_tprs = [data[1] for data in fprtprs]
        rxn_fprtpr[rxn] = fprtprs
        # save it to figure and to file
        #cur_auc  = auc(cur_fprs, cur_tprs, reorder=True)
        #draw_roc_curve(rxn, cur_fprs, cur_tprs, cur_auc):
    # end for ...

    # save it to a file
    save_to_file('../evaluations/rxn_fpr_tprs.txt', rxn_fprtpr)

    print ' .. done '
    #Fin.


"""
    Main (Evaluations)
"""
if __name__ == '__main__':
    evaluations()
    #Fin.