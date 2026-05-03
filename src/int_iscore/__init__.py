"""
int-iScore: Integrated AI Predictions Interface Analysis Scores Toolkit

A comprehensive toolkit for multi-modal protein complex interface analysis.
"""

__version__ = "1.0.0"
__author__ = "Yulin Luo"
__email__ = "luoyl2022@shanghaitech.edu.cn"

try:
    from .core.metrics import (
        procheck_analysis,
        count_clashes,
        calculate_sasa,
        extract_atoms,
        contact_detect,
        dp2_and_cpscore,
        convert_and_soap,
    )
    from .core.parser import (
        parse_config,
    )
    _core_available = True
except ImportError as e:
    _core_available = False
    procheck_analysis = None
    count_clashes = None
    calculate_sasa = None
    extract_atoms = None
    contact_detect = None
    dp2_and_cpscore = None
    convert_and_soap = None
    parse_config = None

try:
    from .utils.analysis import (
        load_csv_data,
        calculate_correlation_weights,
        generate_weights_from_csv,
        analyze_csv_correlations,
        create_int_score_pipeline,
        DEFAULT_WEIGHTS,
        DEFAULT_FEATURE_COLUMNS,
    )
    _utils_available = True
except ImportError as e:
    _utils_available = False
    load_csv_data = None
    calculate_correlation_weights = None
    generate_weights_from_csv = None
    analyze_csv_correlations = None
    create_int_score_pipeline = None
    DEFAULT_WEIGHTS = None
    DEFAULT_FEATURE_COLUMNS = None

__all__ = [
    "procheck_analysis",
    "count_clashes", 
    "calculate_sasa",
    "extract_atoms",
    "contact_detect",
    "dp2_and_cpscore",
    "convert_and_soap",
    "parse_config",
    "load_csv_data",
    "calculate_correlation_weights",
    "generate_weights_from_csv",
    "analyze_csv_correlations",
    "create_int_score_pipeline",
    "DEFAULT_WEIGHTS",
    "DEFAULT_FEATURE_COLUMNS",
    "_core_available",
    "_utils_available",
]
