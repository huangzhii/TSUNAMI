import json
import requests

def enrichr_main(genes_str, gene_set_library):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    # genes_str = '\n'.join([
    #     'PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7',
    #     'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1', 'ANAPC16', 'TMCC1',
    #     'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2', 'PLBD2', 'LARP7', 'TECPR2', 
    #     'ZNF302', 'CUX1', 'MOB2', 'CYTH2', 'SEC22C', 'EIF4E3', 'ROBO2',
    #     'ADAMTS9-AS2', 'CXXC1', 'LINC01314', 'ATF7', 'ATP5F1'
    # ])
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    print("Example Gene List")
    print(data)

    user_list_id = data['userListId']

    ## View added gene list
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'
    response = requests.get(ENRICHR_URL % user_list_id)
    if not response.ok:
        raise Exception('Error getting gene list')
        
    data = json.loads(response.text)
    print("View added gene list")
    print(data)

    ## Get enrichment results
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    # gene_set_library = 'GO_Biological_Process_2017b'
    #'KEGG_2016' 'GO_Molecular_Function_2017b' 'GO_Cellular_Component_2017b' 'GO_Biological_Process_2017b'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    print("Enrichment Results:")
    print(data)
    return(data)
