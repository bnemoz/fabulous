#!/usr/bin/env python3


#############################################################################################################
#                                                                                                           #
#               Imports                                                                                     #
#                                                                                                           #
#############################################################################################################


from dataclasses import dataclass
from enum import Enum
from collections import Counter

import os
import string
import random
import datetime
# import re

import pandas as pd
import numpy as np

import subprocess as sp



#############################################################################################################
#                                                                                                           #
#               Python functions to be used by the API                                                      #
#                                                                                                           #
#############################################################################################################



billing_folder = '/home/serveradmin/Apps/fabulous/billing/'
authentification_folder = '/home/serveradmin/Apps/fabulous/authentication/'


def billing(user, token, app):
    billing_file = os.path.join(billing_folder, f'billing_{user}.csv')
    with open(billing_file, 'a') as f:
        f.write(f"{user}\t{token}\t{app}\t{datetime.datetime.now()}\n")
    return None


def get_bill(user):
    billing_file = os.path.join(billing_folder, f'billing_{user}.csv')
    if os.path.exists(billing_file):
        df = pd.read_csv(billing_file, sep='\t', header=None)
        df.columns=["client ID","Auth token", "App", "Timestamp"]        
        return df
    else:
        return None
    

def authenticate(user, token):
    df = pd.read_csv(os.path.join(authentification_folder, 'credentials.csv'), header=None)
    if token in df[1].values:
        if user == df.loc[df[1] == token, 0].values[0]:
            return True
        else:
            return False
    else:
        return False
