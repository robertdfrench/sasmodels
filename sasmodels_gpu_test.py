def gpu_test():
    import json
    import platform

    import sasmodels
    from sasmodels.model_test import model_tests
    try:
        from sasmodels.kernelcl import environment
        env = environment()
        clinfo = [(ctx.devices[0].platform.vendor,
                      ctx.devices[0].platform.version,
                      ctx.devices[0].vendor,
                      ctx.devices[0].name,
                      ctx.devices[0].version)
                    for ctx in env.context]
    except ImportError:
        clinfo = None

    failures = []
    for test in model_tests():
        try:
            test()
        except Exception:
            failures.append(test.description)

    info = {
        'version':  sasmodels.__version__,
        'platform': platform.uname(),
        'opencl': clinfo,
        'failing tests': failures,
    }
    print(json.dumps(info))

if __name__=="__main__":
    gpu_test()