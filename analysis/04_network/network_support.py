from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/04_network"))
from PCMCI import results_to_long_df
from tigramite import data_processing as pp, plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.gpdc import GPDC
FONT_FAMILY="Arial"
EDGE_COLORS=["#E5BD9B","#CF824A","#9D5024"]
NODE_COLORS=["#EDF3F9","#B4CDDF","#7EA6C6"]
NODE_POSITIONS={2:{"x":[-.65,.65],"y":[0,0]},3:{"x":[-.65,.65,0],"y":[.35,.35,-.50]}}
