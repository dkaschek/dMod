#!/usr/bin/env python3
#
# (c) INCOME Hackathon 2018, Bernried, Daniel^2
#

import sys
import numpy as np
import json

try: 
  import amici.sbml_import
except:
  from amici import sbml_import
  
def symengineMatrixToNumpy(x, astype='float'):
    return np.array(x).reshape(x.shape).astype(astype)

def getModelJSON(sbml_file_name):

    importer = amici.sbml_import.SbmlImporter(sbml_file_name, check_validity=False)


    observables = amici.sbml_import.assignmentRules2observables(importer.sbml,
                                                filter_function=lambda variableId:
    variableId.getId().startswith('observable_') and not
    variableId.getId().endswith('_sigma'))
    importer.processSBML()
    # importer.computeModelEquations()
    
    S = symengineMatrixToNumpy(importer.stoichiometricMatrix)
    dataPy = {
        'S': importer.stoichiometricMatrix.tolist(),
        'v': [str(x) for x in importer.fluxVector],
        'p': importer.parameterIndex,
        'stateNames': symengineMatrixToNumpy(importer.symbols['species']['sym'], astype='str').tolist(),
        'parameterNames': symengineMatrixToNumpy(importer.symbols['parameter']['sym'], astype='str').tolist(),
        'x0': symengineMatrixToNumpy(importer.speciesInitial, astype='str').tolist(),
        "observables": observables
    }
    data = json.dumps(dataPy)

    return data

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s SBML-FILE-NAME [OUTFILE]' % __file__)
        sys.exit(1)

    sbml_file_name = sys.argv[1]
    output = getModelJSON(sbml_file_name)

    if len(sys.argv) > 2:
        outfile = sys.argv[2]
        with open(outfile, "w") as f:
            f.write(output)
    else:
        print(output)
