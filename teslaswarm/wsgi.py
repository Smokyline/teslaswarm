"""
WSGI config for teslaswarm project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/3.0/howto/deployment/wsgi/
"""

import os
from teslaswarm.settings import ACTIVATE_THIS_PATH
from django.core.wsgi import get_wsgi_application
from subprocess import call

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'teslaswarm.settings')

application = get_wsgi_application()

try:
    exec(open(ACTIVATE_THIS_PATH).read())
except Exception as e:
    print('no activate_this.py')