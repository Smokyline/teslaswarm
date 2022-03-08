import os
import sys
import site

# Add the site-packages of the chosen virtualenv to work with
site.addsitedir('home/ivan/anaconda3/envs/tesla_env/lib/python3.7/site-packages')

# Add the app's directory to the PYTHONPATH
sys.path.append('/home/ivan/djangoprojects/teslaswarm')
sys.path.append('/home/ivan/djangoprojects/teslaswarm/teslaswarm')

os.environ['DJANGO_SETTINGS_MODULE'] = 'teslaswarm.settings'

# Activate your virtual env
activate_env=os.path.expanduser("/home/ivan/anaconda3/envs/tesla_env/bin/activate_this.py")
execfile(activate_env, dict(__file__=activate_env))

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
