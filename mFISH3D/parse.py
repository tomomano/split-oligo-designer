from Bio.SeqUtils import MeltingTemp as mt

def oligominer_param_parser(oligominer_param):
    l = oligominer_param['l']
    L = oligominer_param['L']
    gcPercent = oligominer_param['gcPercent']
    GCPercent = oligominer_param['GCPercent']
    tm = oligominer_param['tm']
    TM = oligominer_param['TM']
    X = oligominer_param['X']
    sal = oligominer_param['sal']
    form = oligominer_param['form']
    sp = oligominer_param['sp']
    concA = oligominer_param['concA']
    concB = oligominer_param['concB']
    headerVal = oligominer_param['headerVal']
    bedVal = oligominer_param['bedVal']
    OverlapModeVal = oligominer_param['OverlapModeVal']
    verbocity = oligominer_param['verbocity']
    reportVal = oligominer_param['reportVal']
    debugVal = oligominer_param['debugVal']
    metaVal = oligominer_param['metaVal']
    outNameVal = oligominer_param['outNameVal']

    # Assign concentration variables based on magnitude.
    if concA >= concB:
        conc1 = concA
        conc2 = concB
    else:
        conc1 = concB
        conc2 = concA

    # Retrieve the stack table.
    exec ('nn_table = mt.%s' % oligominer_param['nn_table'])
    return (l, L, gcPercent, GCPercent, nn_table, tm, TM,
                   X, sal, form, sp, conc1, conc2, headerVal, bedVal,
                   OverlapModeVal, verbocity, reportVal, debugVal, metaVal,
                   outNameVal)

