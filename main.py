#use mkdssp to comupet acc
#solution one with error : TypeError: a bytes-like object is required, not 'str'

import tempfile
import subprocess
script_mkdssp = '''\ mkdssp -i 1uaz.pdb -o 1uaz.dssp'''

def run_script(script):
    with tempfile.NamedTemporaryFile() as scriptfile:
        scriptfile.write(script)
        scriptfile.flush()
        subprocess.call(['/bin/bash', scriptfile.name, 'r'])

run_script(script_mkdssp)

#solution two takes only the first argument (mkdssp) when adding more it gives an error
#subprocess.run(['mkdssp' '-i' '1uaz.pdb'])
