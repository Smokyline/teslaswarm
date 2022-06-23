import os
import sys

# Find code directory relative to our directory
sys.path.append(os.getcwd())


def check_matlab_status():
    try:
        import subprocess
        from teslaswarm.settings import BASE_DIR, MATLAB_PATH
        # god pls help us....
        print(MATLAB_PATH)
        cmd = [MATLAB_PATH, 'matlab -nodesktop -nosplash -r "meas_extr_SA(); exit"']
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=BASE_DIR + '/DMA')
        while p.poll() is None:
            l = p.stdout.readline()  # This blocks until it receives a newline.
            print(l)
        print(p.stdout.read())
        print('matlab successfully load')
    except Exception as e:
        print(e)
        print('matlab model load is failed')


if __name__ == '__main__':
    check_matlab_status()
