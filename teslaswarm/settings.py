"""
Django settings for teslaswarm project.

Generated by 'django-admin startproject' using Django 3.0.3.

For more information on this file, see
https://docs.djangoproject.com/en/3.0/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/3.0/ref/settings/
"""

import os
import platform
# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/3.0/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'v@2p7%@z)zns-n-c^1k)1(m*4wiue)gx4za%bvph$l%2d^tdr='

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['aleph.gcras.ru', '127.0.0.1']
print('django.settings ALLOWED_HOSTS', ALLOWED_HOSTS)

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'teslaswarm.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, "templates")],
        #'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'teslaswarm.wsgi.application'


# Database
# https://docs.djangoproject.com/en/3.0/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}


# Password validation
# https://docs.djangoproject.com/en/3.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/3.0/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Moscow'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.0/howto/static-files/




#STATIC_URL = '/static/'
STATIC_URL = '/request/static/'

STATIC_OS_PATH = os.path.join(BASE_DIR, 'static')
print(STATIC_OS_PATH)
#STATIC_OS_PATH = BASE_DIR + '/static'
STATICFILES_DIRS = [
    BASE_DIR + '/static',
    BASE_DIR + "/templates/static",
]
print('django.settings STATICFILES_DIRS', STATICFILES_DIRS)

CHAOS_PATH = os.path.join(BASE_DIR, 'chaos7_model\\data\\CHAOS-7.mat')
#CHAOS_PATH = os.path.join('C:\\Users\\ivan\\YandexDisk\\workspace\\py\\teslaswarm\\chaos7_model\\data\\CHAOS-7.mat')

OBS_DATA_PATH = os.path.join(BASE_DIR, 'static/data/obs_data/SALU/')

ACTIVATE_THIS_PATH = "/home/ivan/anaconda3/bin/activate_this.py"
print('django.settings ACTIVATE_THIS_PATH', ACTIVATE_THIS_PATH)


if platform.system() == 'Linux':
    CDF_LIB_PATH = "/usr/local/cdf/lib/"
elif platform.system() == 'Darwin':
    CDF_LIB_PATH = '/Applications/cdf/cdf38_0-dist/lib/'
elif platform.system() == 'Windows':
    CDF_LIB_PATH = "C:\\CDF_Distribution\\cdf38_0-dist\\lib"
os.environ["CDF_LIB"] = CDF_LIB_PATH
