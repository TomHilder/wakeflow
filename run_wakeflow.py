import sys
from parameters import Parameters, Constants


p = Parameters('config.yaml')
if p.do_sanity_checks() != True:
    print('Exiting')
    sys.exit()
