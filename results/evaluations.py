
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
            if 0 <= idx < 5:
                param_value = each_line[1].strip()
                if isfloat(param_value): param_value = float(param_value)
                ret_data['parameters'][each_line[0]] = param_value
            # control data read
            if idx == 5: 
                data = each_line[1:]
                data = [each_data for each_data in data if each_data]
                ret_data['controls']['ess']     = data
            if idx == 6:
                data = each_line[1:]
                data = [each_data for each_data in data if each_data]
                ret_data['controls']['ess-recon1'] = data
            if idx == 7:
                data = each_line[1:]
                data = [each_data for each_data in data if each_data]
                ret_data['controls']['all-recon1'] = data
            # results data read (from 9th row)
            if idx > 8:
                if idx % 2 == 1:
                    rxn_remove = each_line[0].strip()
                    rxn_remove = rxn_remove.replace('[', '')
                    rxn_remove = rxn_remove.replace(']', '')
                    ret_data['results'][rxn_remove] = {}
                if idx % 2 == 0:
                    data = each_line[1:]
                    data = [each_data for each_data in data if each_data]
                    ret_data['results'][rxn_remove]['ess']     = data
            # end read
    return ret_data

def compute_fprtprs(result):
    ret_data = {
        'parameters' : result['parameters'],
        'results'    : {},
    }
    # control sets
    her2_ess   = set(result['controls']['ess'])
    recon1_all = set(result['controls']['all-recon1'])
    her2_essr  = set(result['controls']['ess-recon1'])
    her2_nessr = recon1_all - her2_essr
    # compute the fpr and tprs
    for each_rxn, each_data in result['results'].iteritems():
        # load the essential list
        cur_ess  = set(each_data['ess'])
        cur_ness = recon1_all - cur_ess
        # compute TP, TN, FP, FN
        true_positive  = len(cur_ess.intersection(her2_essr))
        false_positive = len(cur_ess.intersection(her2_nessr))
        false_negative = len(cur_ness.intersection(her2_essr))
        true_negative  = len(cur_ness.intersection(her2_nessr))
        # compute tpr and fpr
        cur_fpr = float(false_positive) / len(her2_nessr)
        cur_tpr = float(true_positive)  / len(her2_essr)
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

def save_to_file(result):
    cur_params  = result['parameters']
    cur_results = result['results']

    filename = '../evaluations/%s_%s_%s_%s_%s.txt' % \
        (cur_params['method'], cur_params['ess threshold'], \
         cur_params['celline ratio'], cur_params['biomass drop ratio'], \
         cur_params['p value threshold'])
    with open(filename, 'wb') as outputfile:
        outputfile.write('[Reaction removed]: (FP-rate, TP-rate) \n')
        for rxn, data in cur_results.iteritems():
            output_str = '(%03.6f, %03.6f)' % (data[0], data[1])
            outputfile.write('[%40s] : %s \n' % (rxn, output_str))
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

    # draw the ROC curve with computed the AUC in each rxn
    for each_res in res_fptps:
        save_to_file(each_res)
    print ' .. done '
    #Fin.


"""
    Main (Evaluations)
"""
if __name__ == '__main__':
    evaluations()
    #Fin.
