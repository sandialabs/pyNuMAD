# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:40:17 2023

@author: evaande
"""

def get_modal_loss_factors(modal_data):
    loss_factors = dict()
    mat_lf = modal_data['mat_loss_factors']
    modes = modal_data['modes']
    for mk in modes:
        tot_loss = 0.0
        tot_U = 0.0
        md = modes[mk]
        for i, s in enumerate(md['stress']):
            mlf = mat_lf[md['material'][i]]
            e = md['strain'][i]
            vol = md['volume'][i]
            U_el = 0.0
            del_U = 0.0
            for j in range(0,6):
                try:
                    p = s[j]*e[j]
                except:
                    p = float(s[j])*float(e[j])
                U_el = U_el + p
                del_U = del_U + p*mlf[j]
            p = 0.5*vol
            U_el = p*U_el
            del_U = p*del_U
            tot_U = tot_U + U_el
            tot_loss = tot_loss + del_U
        loss_factors[mk] = tot_loss/tot_U
    return loss_factors