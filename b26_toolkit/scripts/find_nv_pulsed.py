from b26_toolkit.scripts.find_nv import FindNV
from b26_toolkit.scripts.set_laser import SetLaser
from b26_toolkit.scripts.galvo_scan.galvo_scan_pulsed import GalvoScanPulsed

class FindNVPulsed(FindNV):
    _SCRIPTS = {'take_image': GalvoScanPulsed, 'set_laser': SetLaser}