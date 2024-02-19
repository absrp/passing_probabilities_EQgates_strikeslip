import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import sklearn
from sklearn import linear_model
from scipy.optimize import curve_fit
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, f1_score, confusion_matrix
from scipy import stats
import warnings
from sklearn.metrics import ConfusionMatrixDisplay
warnings.filterwarnings('ignore')
import matplotlib.colors as mcolors
from scipy.stats import lognorm
from scipy.stats import expon
from scipy.stats import bootstrap
from sklearn.utils import resample


def bootstrap_errors(xfeature, BUbin, classweightb, ax, minx, maxx, length, n_iterations=10000):

    """
    This function calculates errors for the logistic fits based on bootstrapping

    """    

    logistic_reg = []
    np.random.seed(42) # for reproducibility

    # Bootstrap resampling and fitting logistic regression
    for i in range(n_iterations):
        # Resample with replacement
        resampled_xfeature, resampled_BUbin = resample(xfeature, BUbin, replace=True)


    # dealing with small number of unbreached gaps
        if len(resampled_BUbin.unique()) < 2:
            i = i-1 # update i to re-run this iteration
            continue  # Skip the iteration of the loop

        probname = sklearn.linear_model.LogisticRegression(penalty='none', class_weight=classweightb).fit(
            np.atleast_2d(resampled_xfeature).T, resampled_BUbin
        )

        x = np.atleast_2d(np.linspace(minx, maxx, 10000)).T
        logistic_reg.append(probname.predict_proba(x)[:, 1])

    
    percentiles_2_5 = np.percentile(logistic_reg, 2.5, axis=0) # 34
    percentiles_97_5 = np.percentile(logistic_reg, 97.5, axis=0) # 68
    
    xi = np.linspace(minx, maxx, 10000)

    if length==True:
        ax.fill_between(10**xi, percentiles_2_5, percentiles_97_5, color='slategray', alpha=0.2)
    else:
        ax.fill_between(xi, percentiles_2_5, percentiles_97_5, color='slategray', alpha=0.2)
    return logistic_reg


def build_logistic_regression(
        grouped,
        groupid, 
        stress_typeYN, 
        length_or_angle, 
        class_weightb, 
        axesid,
        minx,
        maxx,
        colorline,
        xlabel,
        ptsize
        ):
    
    """
    This function builds logistic regressions for earthquake gates, based on the groups mapped as breached and unbreached.

    """   
    
    EQgate = grouped.get_group(groupid)

    if stress_typeYN  == 'restraining':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('restraining')
        feature = group['Length (m) or angle (deg)']
        BUbin = pd.get_dummies(group['Breached or unbreached'])
        BUbin = BUbin['unbreached']

    elif stress_typeYN  == 'releasing':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('releasing')
        feature = group['Length (m) or angle (deg)']
        BUbin = pd.get_dummies(group['Breached or unbreached'])
        BUbin = BUbin['unbreached']
    
    elif stress_typeYN == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('single')
        feature = group['Length (m) or angle (deg)']
        BUbin = pd.get_dummies(group['Breached or unbreached'])
        BUbin = BUbin['unbreached']

    elif stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        feature = group['Length (m) or angle (deg)']
        BUbin = pd.get_dummies(group['Breached or unbreached'])
        BUbin = BUbin['unbreached']

    else:
        group = EQgate 
        feature = EQgate['Length (m) or angle (deg)']
        group['logfeature'] = np.log10(group['Length (m) or angle (deg)'].astype('float'))
        BUbin = pd.get_dummies(EQgate['Breached or unbreached'])
        BUbin = BUbin['unbreached']
    
    if length_or_angle == 'length':
        group['logfeature'] = np.log10(group['Length (m) or angle (deg)'].astype('float')) 
        xfeature = group['logfeature']
        minx = np.log10(minx)
        maxx = np.log10(maxx)

    elif length_or_angle == 'angle':
        group['logfeature'] = np.log10(group['Length (m) or angle (deg)'])
        xfeature = group['Length (m) or angle (deg)']
        # minx = np.log10(minx)
        # maxx = np.log10(maxx)
    else: 
        KeyError("Feature must include a length or an angle")


    palette = {'breached': 'teal', 'unbreached': 'darkorange'}
    
    if max(group['Length (m) or angle (deg)'])>90:
            sns.swarmplot(
            data=group,
            #x = 'logfeature',
            x='Length (m) or angle (deg)',
            y='Breached or unbreached',
            ax=axesid,size=ptsize,
            hue="Breached or unbreached",
            palette=palette,alpha=0.7
        )
    else:  
        sns.swarmplot(
            data=group,
            x='Length (m) or angle (deg)',
            y='Breached or unbreached',
            ax=axesid,size=ptsize,
            hue="Breached or unbreached",
            palette=palette,alpha=0.7
        )

    if max(group['Length (m) or angle (deg)'])>90:
        bootstrap_errors(xfeature, BUbin, class_weightb, axesid, minx, maxx, True)
    else: 
        bootstrap_errors(xfeature, BUbin, class_weightb, axesid, minx, maxx, False)

    probname = sklearn.linear_model.LogisticRegression(penalty='none',class_weight=class_weightb).fit(np.atleast_2d(xfeature).T,BUbin)

    # tests
    acc = accuracy_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    pre = precision_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    f1 = f1_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    roc =  roc_auc_score(BUbin, probname.predict_proba(np.atleast_2d(xfeature).T)[:,1])
    confusion_matrixi = confusion_matrix(BUbin, probname.predict(np.atleast_2d(xfeature).T))

    x = np.atleast_2d(np.linspace(minx, maxx, 10000)).T

    if max(group['Length (m) or angle (deg)'])>90:
        axesid.plot(10**x,probname.predict_proba(x)[:,1],color = colorline)
        axesid.text(10**x[-10], -0.1, f'ROC={roc:.2f}', ha='right', va='top',fontsize=14)
        axesid.set_xscale('log')
    else:
        axesid.text(x[-10], -0.1, f'ROC={roc:.2f}', ha='right', va='top',fontsize=14)
        # axesid.text(x[-10], 0.1, r'$P_{50}$=' + str(int(x_at_threshold)) + '$^{\circ}$', ha='right', va='top', fontsize=14)
        axesid.plot(x,probname.predict_proba(x)[:,1],color = colorline)

    axesid.set_ylabel('Passing probability')
    axesid.set_xlabel(xlabel)
    axesid.set_yticklabels(["Breached", "Unbreached"],rotation=90,va='center')
    axesid.get_legend().remove()

    return probname, acc, pre, f1, roc, confusion_matrixi, BUbin, xfeature 

def build_regression_double_bend_length(grouped,groupid,stress_typeYN,feature_type,axesid,minx,maxx,xlabel,ptsize,colorline='slategrey',class_weightb=None):

    """
    This function calculates a logistic regression for a set of earthquake gate double bend length metrics from surface ruptures and plots the results.

    """    
    EQgate = grouped.get_group(groupid)

    if stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        feature = group[feature_type]
        BUbin = pd.get_dummies(group['Breached or unbreached'])
        BUbin = BUbin['breached']
        #group['logfeature'] = np.log10(group['Length (m) or angle (deg)'].astype('float')) 
        xfeature = group['Length (m) or angle (deg)'].astype('float')
        minx = minx
        maxx = maxx

    palette = {'breached': 'teal', 'unbreached': 'darkorange'}
    if max(xfeature)>90:
        sns.swarmplot(data=group,x=feature_type,y='Breached or unbreached',ax=axesid,size=ptsize,hue="Breached or unbreached",palette=palette,alpha=0.7)
    else:
        sns.swarmplot(data=group,x=feature_type,y='Breached or unbreached',ax=axesid,size=ptsize,hue="Breached or unbreached",palette=palette,alpha=0.7)

    axesid.set_xscale('log')
    axesid.set_xlabel(xlabel)
    axesid.set_yticklabels(["Breached", "Unbreached"],rotation=90,va='center')
    axesid.get_legend().remove()

    if max(group['Length (m) or angle (deg)'])>90:
        bootstrap_errors(xfeature, BUbin, class_weightb, axesid, minx, maxx, True)
    else: 
        bootstrap_errors(xfeature, BUbin, class_weightb, axesid, minx, maxx, False)

    probname = sklearn.linear_model.LogisticRegression(penalty='none',class_weight=class_weightb).fit(np.atleast_2d(xfeature).T,BUbin)

    # tests
    acc = accuracy_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    pre = precision_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    f1 = f1_score(BUbin, probname.predict(np.atleast_2d(xfeature).T))
    roc =  roc_auc_score(BUbin, probname.predict_proba(np.atleast_2d(xfeature).T)[:,1])
    confusion_matrixi = confusion_matrix(BUbin, probname.predict(np.atleast_2d(xfeature).T))

    x = np.atleast_2d(np.linspace(minx, maxx, 10000)).T

    if max(group['Length (m) or angle (deg)'])>90:
        axesid.plot(x,probname.predict_proba(x)[:,1],color = colorline)
        axesid.text(x[-10], -0.1, f'ROC={roc:.2f}', ha='right', va='top',fontsize=14)
        axesid.set_xscale('log')
    else:
        axesid.text(x[-10], -0.1, f'ROC={roc:.2f}', ha='right', va='top',fontsize=14)
        # axesid.text(x[-10], 0.1, r'$P_{50}$=' + str(int(x_at_threshold)) + '$^{\circ}$', ha='right', va='top', fontsize=14)
        axesid.plot(x,probname.predict_proba(x)[:,1],color = colorline)

    axesid.set_ylabel('Passing probability')
    axesid.set_xlabel(xlabel)
    axesid.set_yticklabels(["Breached", "Unbreached"],rotation=90,va='center')

    return BUbin, xfeature 

# make CDFs
def build_cdf(grouped,groupid, stress_typeYN, length_or_angle, colorB,colorU, axesid,xlabel,labelB, labelU):

    """
    This function generates CDFs for the geometry measured for each earthquake gate, separated into breached and unbreached features and restraining and releasing when available. 

    """
    EQgate = grouped.get_group(groupid)

    if stress_typeYN  == 'restraining':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('restraining')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN  == 'releasing':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('releasing')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('single')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif groupid == 'strand':
        grouped_BU = EQgate.groupby(EQgate["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('breached')

    else: 
        group = EQgate
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
        
    if length_or_angle == 'length':
        xvals_U = np.log10(xvals_U['Length (m) or angle (deg)'])
        xvals_B = np.log10(xvals_B['Length (m) or angle (deg)'])
        #xvals_U = np.log10(xvals_U['Spacing double bend (m)'])
        #xvals_B = np.log10(xvals_B['Spacing double bend (m)'])
        

    elif length_or_angle == 'angle':
        xvals_B = xvals_B['Length (m) or angle (deg)']
        xvals_U = xvals_U['Length (m) or angle (deg)']

    else: 
        KeyError("Feature must include a length or an angle")

    # sort the data:
    sortB = np.sort(xvals_B)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(xvals_B)) / (len(xvals_B) - 1) 
    sns.ecdfplot(sortB,c=colorB,label=labelB,ax=axesid) 
    axesid.set_xlabel(xlabel)
    sortU = np.sort(xvals_U)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(xvals_U)) / (len(xvals_U) - 1) 
    #axesid.plot(sortU,p,c=colorU,label=labelU) 
    sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
    axesid.set_xlabel(xlabel)
    axesid.set_ylabel('Proportion')
    axesid.grid(color='lightgray', linewidth=0.5, alpha=0.5)
    axesid.legend()


def build_cdf_lognorm(grouped, groupid, stress_typeYN, length_or_angle, colorB, colorU, axesid,xlabel,labelB, labelU):

    """
    This function generates CDFs for the geometry measured for each earthquake gate, separated into breached and unbreached features and restraining and releasing when available. We then fit log normal CDFs to the ECDF.

    """
    EQgate = grouped.get_group(groupid)

    if stress_typeYN  == 'restraining':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('restraining')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN  == 'releasing':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('releasing')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('single')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif groupid == 'strand':
        grouped_BU = EQgate.groupby(EQgate["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('breached')

    else: 
        group = EQgate
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
        
    if length_or_angle == 'length':
        xvals_U = xvals_U['Length (m) or angle (deg)']
        xvals_B = xvals_B['Length (m) or angle (deg)']
        #xvals_U = np.log10(xvals_U['Spacing double bend (m)'])
        #xvals_B = np.log10(xvals_B['Spacing double bend (m)'])
        

    elif length_or_angle == 'angle':
        xvals_B = xvals_B['Length (m) or angle (deg)']
        xvals_U = xvals_U['Length (m) or angle (deg)']

    else: 
        KeyError("Feature must include a length or an angle")


    sortB = np.sort(xvals_B)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(xvals_B)) / (len(xvals_B) - 1) 
    sns.ecdfplot(sortB,c=colorB,label=labelB,ax=axesid) 
    axesid.set_xlabel(xlabel)

    # log normal fit 
    yvals = np.arange(len(sortB)) / float(len(sortB)-1)
    shape, loc, scale = lognorm.fit(sortB, floc=0, f0=1-yvals[-1])
    xvals = np.linspace(min(sortB), max(sortB), 100)
    cdf_fitted = lognorm.cdf(xvals, shape, loc, scale)
    axesid.plot(xvals, cdf_fitted, color=colorB,linestyle=':')  

    sortU = np.sort(xvals_U)
    yvals = np.arange(len(sortU)) / float(len(sortU)-1)
    shape, loc, scale = lognorm.fit(sortU, floc=0, f0=1-yvals[-1])
    xvals = np.linspace(min(sortU), max(sortU), 100)
    cdf_fitted = lognorm.cdf(xvals, shape, loc, scale)

    if max(xvals)>90:
        sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
        axesid.plot(xvals, cdf_fitted, color=colorU,linestyle=':') 
        axesid.set_xscale('log')
    else: 
        sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
        axesid.plot(xvals, cdf_fitted, color=colorU,linestyle=':') 
    
    # exponential fit
    # yvals = np.arange(len(sortB)) / float(len(sortB)-1)
    # loc, scale = expon.fit(sortB, floc=0)
    # xvals = np.linspace(min(sortB), max(sortB), 100)
    # cdf_fitted = expon.cdf(xvals, loc, scale)
    # axesid.plot(xvals, cdf_fitted, color=colorB,linestyle='--')  

    # sortU = np.sort(xvals_U)
    # yvals = np.arange(len(sortU)) / float(len(sortU)-1)
    # loc, scale = expon.fit(sortU, floc=0)
    # xvals = np.linspace(min(sortU), max(sortU), 100)
    # cdf_fitted = expon.cdf(xvals, loc, scale)

    # if max(xvals)>90:
    #     sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
    #     axesid.plot(xvals, cdf_fitted, color=colorU,linestyle='--') 
    #     axesid.set_xscale('log')
    # else: 
    #     sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
    #     axesid.plot(xvals, cdf_fitted, color=colorU,linestyle='--') 

    axesid.set_xlabel(xlabel)
    axesid.set_ylabel('Proportion')
    if groupid == 'stepover':
        axesid.legend(fontsize=8)
    else:
        axesid.legend(fontsize=10)

    axesid.grid(color='lightgray', linewidth=0.5, alpha=0.5)


# make PDFs
def build_pdf(grouped,groupid, stress_typeYN, length_or_angle, colorB,colorU, axesid,xlabel,labelB, labelU):

    """
    This function calculates PFDs for each earthquake gate geometry.
    """
    EQgate = grouped.get_group(groupid)

    if stress_typeYN  == 'restraining':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('restraining')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN  == 'releasing':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('releasing')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('single')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif groupid == 'strand':
        grouped_BU = EQgate.groupby(EQgate["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('breached')

    else: 
        group = EQgate
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
        
    if length_or_angle == 'length':
        xvals_U = np.log10(xvals_U['Length (m) or angle (deg)'])
        xvals_B = np.log10(xvals_B['Length (m) or angle (deg)'])
        

    elif length_or_angle == 'angle':
        xvals_B = xvals_B['Length (m) or angle (deg)']
        xvals_U = xvals_U['Length (m) or angle (deg)']

    else: 
        KeyError("Feature must include a length or an angle")

    sns.kdeplot(xvals_U,c=colorU,label=labelU,ax=axesid)
    sns.kdeplot(xvals_B,c=colorB,label=labelB,ax=axesid)
    
    axesid.set_xlabel(xlabel)
    axesid.set_ylabel('Proportion')
    axesid.legend()


# make CDFs bend lengths
def build_cdf_bend_lengths(grouped,groupid, stress_typeYN, featuretype, colorB,colorU, axesid,xlabel,labelB, labelU):

    """

    This function generates CDFs for the bend length and proxy step-over length for double bends, separated into breached and unbreached groups

    """
    EQgate = grouped.get_group(groupid)

    if stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
    
        xvals_U = xvals_U[featuretype]
        xvals_B = xvals_B[featuretype]


    else: 
        KeyError("This function only works for double bends")

    # sort the data:
    sortB = np.sort(xvals_B)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(xvals_B)) / (len(xvals_B) - 1) 
    sns.ecdfplot(sortB,c=colorB,label=labelB,ax=axesid) 
    axesid.set_xlabel(xlabel)
    sortU = np.sort(xvals_U)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(xvals_U)) / (len(xvals_U) - 1) 
    #axesid.plot(sortU,p,c=colorU,label=labelU) 
    sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
    axesid.set_xlabel(xlabel)
    axesid.set_xscale('log')
    axesid.set_ylabel('Proportion')
    axesid.legend()
    axesid.legend(fontsize=10)


    # log normal fit breached
    yvals = np.arange(len(sortB)) / float(len(sortB)-1)
    shape, loc, scale = lognorm.fit(sortB, floc=0, f0=1-yvals[-1])
    xvals = np.linspace(min(sortB), max(sortB), 100)
    cdf_fitted = lognorm.cdf(xvals, shape, loc, scale)
    axesid.plot(xvals, cdf_fitted, color=colorB,linestyle=':')  

    sortU = np.sort(xvals_U)
    yvals = np.arange(len(sortU)) / float(len(sortU)-1)
    shape, loc, scale = lognorm.fit(sortU, floc=0, f0=1-yvals[-1])
    xvals = np.linspace(min(sortU), max(sortU), 100)
    cdf_fitted = lognorm.cdf(xvals, shape, loc, scale)

    if max(xvals)>90:
        sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
        axesid.plot(xvals, cdf_fitted, color=colorU,linestyle=':') 
        axesid.set_xscale('log')
    else: 
        sns.ecdfplot(sortU,c=colorU,label=labelU,ax=axesid)
        axesid.plot(xvals, cdf_fitted, color=colorU,linestyle=':') 
    axesid.grid(color='lightgray', linewidth=0.5, alpha=0.5)

def kstest_variables(grouped,groupid, stress_typeYN, length_or_angle):

    """
    This function runs a ks test through two populations to compare whether they are drawn from the same distribution.
    """

    EQgate = grouped.get_group(groupid)
    
    if stress_typeYN  == 'restraining':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('restraining')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN  == 'releasing':
        grouped_stress = EQgate.groupby(EQgate["Type (releasing or restraining)"])
        group = grouped_stress.get_group('releasing')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('single')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')

    elif stress_typeYN == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        group = grouped_stress.get_group('double')
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
    
    elif stress_typeYN == 'releasing_restraining_breached':
        grouped_stress = EQgate.groupby(EQgate["Breached or unbreached"])
        group = grouped_stress.get_group('breached')
        grouped_BU = group.groupby(group["Type (releasing or restraining)"])
        xvals_B = grouped_BU.get_group('releasing')
        xvals_U = grouped_BU.get_group('restraining')

    elif stress_typeYN == 'releasing_restraining_unbreached':
        grouped_stress = EQgate.groupby(EQgate["Breached or unbreached"])
        group = grouped_stress.get_group('unbreached')
        grouped_BU = group.groupby(group["Type (releasing or restraining)"])
        xvals_B = grouped_BU.get_group('releasing')
        xvals_U = grouped_BU.get_group('restraining')

    elif groupid == 'strand':
        grouped_BU = EQgate.groupby(EQgate["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('breached')

    else: 
        group = EQgate
        grouped_BU = group.groupby(group["Breached or unbreached"])
        xvals_B = grouped_BU.get_group('breached')
        xvals_U = grouped_BU.get_group('unbreached')
        
    if length_or_angle == 'length':
        xvals_U = np.log10(xvals_U['Length (m) or angle (deg)'])
        xvals_B = np.log10(xvals_B['Length (m) or angle (deg)'])
 #       xvals_U = np.log10(xvals_U['Spacing double bend (m)'])
 #       xvals_B = np.log10(xvals_B['Spacing double bend (m)'])

    elif length_or_angle == 'angle':
        xvals_B = xvals_B['Length (m) or angle (deg)']
        xvals_U = xvals_U['Length (m) or angle (deg)']

    else: 
        KeyError("Feature must include a length or an angle")
    
    return  stats.kstest(xvals_B, xvals_U)

# # Define the power-law function
def power_law(x, a, b):
     return b*np.log10(x)+a

def gate_distribution_along_strike(grouped,groupid,type,length_or_angle,axesid,ylabel,palette):

    """
    This function plots the distribution of earthquake gates of a given type along the surface rupture
    
    """
        
    EQgate = grouped.get_group(groupid)

    if type == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        single = grouped_stress.get_group('single')
        feature = single['Length (m) or angle (deg)']
        normalized_loc = single['Normalized location']

    elif type == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        double = grouped_stress.get_group('double')
        feature = double['Length (m) or angle (deg)']
        normalized_loc = double['Normalized location']

    else:
        feature = EQgate['Length (m) or angle (deg)']
        normalized_loc = EQgate['Normalized location']
    
    if length_or_angle == 'length':
        yfeature = feature

    elif length_or_angle == 'angle':
        yfeature = feature

    else: 
        KeyError("Feature must include a length or an angle")

    # now we plot the data and fit the model:    
    sns.scatterplot(data=EQgate,x=normalized_loc,y=yfeature,hue=EQgate['Breached or unbreached'],palette=palette, edgecolor='none',alpha=0.6, ax=axesid,legend='')
    axesid.set_ylabel(ylabel)
    if max(yfeature)>90:
        axesid.set_yscale('log')
    axesid.set_xlabel('Normalized distance along the rupture')

def gate_distribution_along_strike_histogram(grouped,groupid,type,length_or_angle,axesid,ylabel,palette):

    """
    This function plots the distribution of earthquake gates of a given type along the surface rupture
    
    """
        
    EQgate = grouped.get_group(groupid)

    if type == 'single':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        single = grouped_stress.get_group('single')
        feature = single['Length (m) or angle (deg)']
        normalized_loc = single['Normalized location']

    elif type == 'double':
        grouped_stress = EQgate.groupby(EQgate["Type (single or double)"])
        double = grouped_stress.get_group('double')
        feature = double['Length (m) or angle (deg)']
        normalized_loc = double['Normalized location']

    else:
        feature = EQgate['Length (m) or angle (deg)']
        normalized_loc = EQgate['Normalized location']
    
    if length_or_angle == 'length':
        yfeature = feature

    elif length_or_angle == 'angle':
        yfeature = feature

    else: 
        KeyError("Feature must include a length or an angle")

 
    sns.histplot(data=EQgate,x=normalized_loc,hue=EQgate['Breached or unbreached'],palette=palette, edgecolor='none',alpha=0.6, ax=axesid,legend='')
    axesid.set_title(ylabel)
    axesid.set_ylabel('Frequency')
    axesid.set_xlabel('Normalized distance along the rupture')

def calculate_center(numbers):
    center_values = []
    for i in range(len(numbers) - 1):
        center = (numbers[i] + numbers[i + 1]) / 2.0
        center_values.append(center)
    return center_values