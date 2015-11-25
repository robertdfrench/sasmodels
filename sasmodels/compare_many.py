#!/usr/bin/env python
import sys
import traceback

import numpy as np

from . import core
from .kernelcl import environment
from .compare import (MODELS, randomize_model, suppress_pd, eval_sasview,
                      eval_opencl, eval_ctypes, make_data, get_demo_pars,
                      columnize, constrain_pars)

def get_stats(target, value, index):
    resid = abs(value-target)[index]
    relerr = resid/target[index]
    srel = np.argsort(relerr)
    #p90 = int(len(relerr)*0.90)
    p95 = int(len(relerr)*0.95)
    maxrel = np.max(relerr)
    rel95 = relerr[srel[p95]]
    maxabs = np.max(resid[srel[p95:]])
    maxval = np.max(value[srel[p95:]])
    return maxrel,rel95,maxabs,maxval

def print_column_headers(pars, parts):
    stats = list('Max rel err|95% rel err|Max abs err above 90% rel|Max value above 90% rel'.split('|'))
    groups = ['']
    for p in parts:
        groups.append(p)
        groups.extend(['']*(len(stats)-1))
    groups.append("Parameters")
    columns = ['Seed'] + stats*len(parts) +  list(sorted(pars.keys()))
    print(','.join('"%s"'%c for c in groups))
    print(','.join('"%s"'%c for c in columns))

def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5):
    model_definition = core.load_model_definition(name)
    pars = get_demo_pars(model_definition)
    header = '\n"Model","%s","Count","%d"'%(name, N)
    if not mono: header += ',"Cutoff",%g'%(cutoff,)
    print(header)

    def trymodel(fn, *args, **kw):
        try:
            result, _ = fn(model_definition, pars_i, data, *args, **kw)
        except KeyboardInterrupt:
            raise
        except:
            print >>sys.stderr, traceback.format_exc()
            print >>sys.stderr, "when comparing",name,"for seed",seed
            if hasattr(data, 'qx_data'):
                result = np.NaN*data.data
            else:
                result = np.NaN*data.x
        return result

    num_good = 0
    first = True
    max_diff = 0
    for k in range(N):
        print >>sys.stderr, name, k
        pars_i, seed = randomize_model(pars)
        constrain_pars(model_definition, pars_i)
        if mono: suppress_pd(pars_i)

        good = True
        labels = []
        columns = []
        if 1:
            sasview_value = trymodel(eval_sasview)
        if 0:
            gpu_single_value = trymodel(eval_opencl, dtype='single', cutoff=cutoff)
            stats = get_stats(sasview_value, gpu_single_value, index)
            columns.extend(stats)
            labels.append('GPU single')
            max_diff = max(max_diff, stats[0])
            good = good and (stats[0] < 5e-5)
        if 0 and environment().has_double:
            gpu_double_value = trymodel(eval_opencl, dtype='double', cutoff=cutoff)
            stats = get_stats(sasview_value, gpu_double_value, index)
            columns.extend(stats)
            labels.append('GPU double')
            max_diff = max(max_diff, stats[0])
            good = good and (stats[0] < 1e-12)
        if 1:
            cpu_double_value = trymodel(eval_ctypes, dtype='double', cutoff=cutoff)
            stats = get_stats(sasview_value, cpu_double_value, index)
            columns.extend(stats)
            labels.append('CPU double')
            max_diff = max(max_diff, stats[0])
            good = good and (stats[0] < 1e-12)
        if 0:
            stats = get_stats(cpu_double_value, gpu_single_value, index)
            columns.extend(stats)
            labels.append('single/double')
            max_diff = max(max_diff, stats[0])
            good = good and (stats[0] < 5e-5)

        columns += [v for _,v in sorted(pars_i.items())]
        if first:
            print_column_headers(pars_i, labels)
            first = False
        if good:
            num_good += 1
        else:
            print(("%d,"%seed)+','.join("%g"%v for v in columns))
    print '"good","%d/%d","max diff",%g'%(num_good, N, max_diff)


def print_usage():
    print "usage: compare_many.py MODEL COUNT (1dNQ|2dNQ) (CUTOFF|mono)"


def print_models():
    print(columnize(MODELS, indent="  "))


def print_help():
    print_usage()
    print("""\

MODEL is the model name of the model or "all" for all the models
in alphabetical order.

COUNT is the number of randomly generated parameter sets to try. A value
of "10000" is a reasonable check for monodisperse models, or "100" for
polydisperse models.   For a quick check, use "100" and "5" respectively.

NQ is the number of Q values to calculate.  If it starts with "1d", then
it is a 1-dimensional problem, with log spaced Q points from 1e-3 to 1.0.
If it starts with "2d" then it is a 2-dimensional problem, with linearly
spaced points Q points from -1.0 to 1.0 in each dimension. The usual
values are "1d100" for 1-D and "2d32" for 2-D.

CUTOFF is the cutoff value to use for the polydisperse distribution. Weights
below the cutoff will be ignored.  Use "mono" for monodisperse models.  The
choice of polydisperse parameters, and the number of points in the distribution
is set in compare.py defaults for each model.

Available models:
""")
    print_models()

def main():
    if len(sys.argv) == 1:
        print_help()
        sys.exit(1)

    model = sys.argv[1]
    if not (model in MODELS) and (model != "all"):
        print 'Bad model %s.  Use "all" or one of:'
        print_models()
        sys.exit(1)
    try:
        count = int(sys.argv[2])
        is2D = sys.argv[3].startswith('2d')
        assert sys.argv[3][1] == 'd'
        Nq = int(sys.argv[3][2:])
        mono = sys.argv[4] == 'mono'
        cutoff = float(sys.argv[4]) if not mono else 0
    except:
        print_usage()
        sys.exit(1)

    data, index = make_data(qmax=1.0, is2D=is2D, Nq=Nq)
    model_list = [model] if model != "all" else MODELS
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono, cutoff=cutoff)

if __name__ == "__main__":
    main()