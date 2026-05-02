from setuptools import setup

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError as exc:
    raise RuntimeError("pybind11 is required to build the sc_backend extension") from exc


ext_modules = [
    Pybind11Extension(
        "int_iscore.utils.sc_backend",
        ["src/int_iscore/utils/sc_backend.cpp"],
        cxx_std=17,
    ),
]


setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
