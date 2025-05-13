def get_worms_from_scientific_name(tax_df, ordered_rank_columns, queue,full_tax_column="taxonomy", like=False, marine_only=False,full_tax_vI = False):

    """
    Using the WoRMS REST API and the pyworm library, retrieve WoRMS ID and taxon ID given a scientific name.

    Dependencies:
        import pyworms, multiprocess, pandas

    Usage:
        get_worms_from_scientific_name(stax_df, ordered_rank_columns, queue,full_tax_column, like, marine_only,full_tax_vI)

    Inputs:
        The scientific name of interest as a string, e.g. 'Dosidicus gigas'
        Optionally, the verbose flag to print species names that aren't matched

    Outputs:
        1. scientificName: WoRMS specified scientific name that matched to sci_name
        2. scientificNameID: WoRMS specified ID
        3. taxonID: WoRMS specified taxon ID

    Katherine Silliman
    2023-10-03
    Python 3.11
    """
# imports
    import pyworms
    import multiprocess as mp
    import pandas as pd

    matches = []
    w_ranks = ['kingdom','phylum','class','order','family','genus']
    for index, row in tax_df.iterrows():
        full_tax = row[full_tax_column]
        row_data = {'full_tax':full_tax,'verbatimIdentification': full_tax}
        if full_tax_vI:
            row_data = {'full_tax':full_tax,'verbatimIdentification': full_tax}
        else:   
            row_data = {'full_tax':full_tax,'verbatimIdentification': 'Null'}
        for i in ordered_rank_columns:
            rank = i
            old_name = row[i]
            if pd.isna(old_name):
                continue 
            else:
                row_data.update({'old_taxonRank': rank, 'old name': old_name})
                if row_data['verbatimIdentification'] == 'Null':
                    row_data['verbatimIdentification'] = old_name
                s_match = pyworms.aphiaRecordsByName(old_name,like=like,marine_only=marine_only)
                #time.sleep(1)
                if s_match == None:
                    row_data['scientificName'] = "No match"
                    row_data['scientificNameID'] = "None"
                    print(old_name+": No match, "+rank)
                    continue
                elif len(s_match) > 1:
                    mult = []
                    for m in s_match:
                        if m['status'] == 'accepted':
                            mult = mult + [m]
                    if len(mult) > 1:
                        row_data['scientificName'] = "Multiple matches"
                        row_data['scientificNameID'] = "None"
                        #print(old_name+": Multiple matches, "+rank+" ")
                    elif len(mult) < 1:
                        row_data['scientificName'] = "Multiple unaccepted matches"
                        row_data['scientificNameID'] = "None"
                        #print(old_name+": Multiple unaccepted matches, "+rank+" ")
                    elif len(mult) == 1:
                        row_data['scientificName'] = mult[0]['scientificname']
                        row_data['scientificNameID'] = mult[0]['lsid']
                        row_data.update(dict(zip(w_ranks, [mult[0].get(key) for key in w_ranks])))
                        row_data.update({'taxonRank': mult[0]['rank']})
                        break
                elif len(s_match) == 1:
                    if s_match[0]['status'] == 'accepted':
                        row_data['scientificName'] = s_match[0]['scientificname']
                        row_data['scientificNameID'] = s_match[0]['lsid']
                        row_data.update(dict(zip(w_ranks, [s_match[0].get(key) for key in w_ranks])))
                        row_data.update({'taxonRank': s_match[0]['rank']})
                        break
                    elif s_match[0]['status'] == 'unaccepted':
                        valid_name = s_match[0]['valid_name']
                        if valid_name != None:
                            v_match = pyworms.aphiaRecordsByName(valid_name,like=like,marine_only=marine_only)
                            row_data['scientificName'] = v_match[0]['scientificname']
                            row_data['scientificNameID'] = v_match[0]['lsid']
                            row_data.update(dict(zip(w_ranks, [v_match[0].get(key) for key in w_ranks])))
                            row_data.update({'taxonRank': v_match[0]['rank']})
                            print(old_name+": Unaccepted, using "+valid_name+", "+rank+" ")
                        #else:
                            #print(old_name+": Unaccepted, no valid name, "+rank+" ")
        matches += [row_data]
    matches = pd.DataFrame.from_dict(matches)
    queue.put(matches)

def get_worms_from_scientific_name_parallel(tax_df, ordered_rank_columns, full_tax_column="taxonomy",like=False, marine_only=False,full_tax_vI = False,n_proc=0):
    import multiprocess as mp
    import pandas as pd
    import os
    
    # Use spawn method on Windows to avoid issues with forking
    if os.name == 'nt':  # Windows
        mp.set_start_method('spawn', force=True)
    
    if n_proc == 0:
        # create as many processes as there are CPUs on your machine
        num_processes = mp.cpu_count()
    else:
        num_processes = n_proc
    
    # calculate the chunk size as an integer
    chunk_size = max(1, round(tax_df.shape[0]/num_processes))
    
    # Create chunks
    chunks = [tax_df.iloc[tax_df.index[i:i + chunk_size]] for i in range(0, tax_df.shape[0], chunk_size)]
    
    # Create a manager to handle the queue
    manager = mp.Manager()
    queue = manager.Queue()
    
    # Create and start processes
    procs = []
    for chunk in chunks:
        proc = mp.Process(
            target=get_worms_from_scientific_name,
            args=(chunk, ordered_rank_columns, queue, full_tax_column, like, marine_only, full_tax_vI)
        )
        procs.append(proc)
        proc.start()
    
    # Collect results from the queue
    results = []
    for _ in range(len(procs)):
        try:
            results.append(queue.get(timeout=600))  # 10-minute timeout
        except Exception as e:
            print(f"Error getting results from queue: {e}")
    
    # Wait for all processes to finish
    for proc in procs:
        proc.join(timeout=10)  # 10-second timeout for joining
        if proc.is_alive():
            print(f"Warning: Process {proc.pid} did not terminate properly, forcing termination")
            proc.terminate()
    
    # Combine results
    new_df = pd.DataFrame()
    for result in results:
        new_df = pd.concat([new_df, result])
    
    return new_df

def get_worms_from_aphiaid_or_name(tax_df, worms_dict, ordered_rank_columns, queue, full_tax_column="taxonomy", like=False, marine_only=False, full_tax_vI=False):
    
    import pyworms
    import multiprocess as mp
    import pandas as pd

    matches = []
    w_ranks = ['kingdom','phylum','class','order','family','genus']
    for index, row in tax_df.iterrows():
        full_tax = row[full_tax_column]
        if row['species'] in worms_dict.keys():
            if full_tax_vI:
                row_data = {'full_tax':full_tax,'verbatimIdentification': full_tax}
            else:
                row_data = {'full_tax':full_tax,'verbatimIdentification': row['species']}
            aid = worms_dict[row['species']]
            record = pyworms.aphiaRecordByAphiaID(aid)
            row_data.update({'taxonRank': 'species', 'old name': 'aphiaID'})
            if record['status'] == 'accepted':
                row_data['scientificName'] = record['scientificname']
                row_data['scientificNameID'] = record['lsid']
                row_data.update(dict(zip(w_ranks, [record.get(key) for key in w_ranks])))
                row_data.update({'taxonRank': record['rank']})
            elif record['status'] == 'unaccepted':
                valid_name = record['valid_name']
                if valid_name != None:
                    v_match = pyworms.aphiaRecordsByName(valid_name,like=like,marine_only=marine_only)
                    row_data['scientificName'] = v_match[0]['scientificname']
                    row_data['scientificNameID'] = v_match[0]['lsid']
                    row_data.update(dict(zip(w_ranks, [v_match[0].get(key) for key in w_ranks])))
                    row_data.update({'taxonRank': v_match[0]['rank']})
                    print(aid+": Unaccepted, using "+valid_name)
                else:
                    print(aid+": Unaccepted, no valid name ")
        else:
            if full_tax_vI:
                row_data = {'full_tax':full_tax,'verbatimIdentification': full_tax}
            else:   
                row_data = {'full_tax':full_tax,'verbatimIdentification': 'Null'}
            for i in ordered_rank_columns:
                rank = i
                old_name = row[i]
                if pd.isna(old_name):
                    continue 
                else:
                    row_data.update({'old_taxonRank': rank, 'old name': old_name})
                    if row_data['verbatimIdentification'] == 'Null':
                        row_data['verbatimIdentification'] = old_name
                    s_match = pyworms.aphiaRecordsByName(old_name,like=like,marine_only=marine_only)
                    #time.sleep(1)
                    if s_match == None:
                        row_data['scientificName'] = "No match"
                        row_data['scientificNameID'] = "None"
                        print(old_name+": No match, "+rank)
                        continue
                    elif len(s_match) > 1:
                        mult = []
                        for m in s_match:
                            if m['status'] == 'accepted':
                                mult = mult + [m]
                        if len(mult) > 1:
                            row_data['scientificName'] = "Multiple matches"
                            row_data['scientificNameID'] = "None"
                            #print(old_name+": Multiple matches, "+rank+" ")
                        elif len(mult) < 1:
                            row_data['scientificName'] = "Multiple unaccepted matches"
                            row_data['scientificNameID'] = "None"
                            #print(old_name+": Multiple unaccepted matches, "+rank+" ")
                        elif len(mult) == 1:
                            row_data['scientificName'] = mult[0]['scientificname']
                            row_data['scientificNameID'] = mult[0]['lsid']
                            row_data.update(dict(zip(w_ranks, [mult[0].get(key) for key in w_ranks])))
                            row_data.update({'taxonRank': mult[0]['rank']})
                            break
                    elif len(s_match) == 1:
                        if s_match[0]['status'] == 'accepted':
                            row_data['scientificName'] = s_match[0]['scientificname']
                            row_data['scientificNameID'] = s_match[0]['lsid']
                            row_data.update(dict(zip(w_ranks, [s_match[0].get(key) for key in w_ranks])))
                            row_data.update({'taxonRank': s_match[0]['rank']})
                            break
                        elif s_match[0]['status'] == 'unaccepted':
                            valid_name = s_match[0]['valid_name']
                            if valid_name != None:
                                v_match = pyworms.aphiaRecordsByName(valid_name,like=like,marine_only=marine_only)
                                row_data['scientificName'] = v_match[0]['scientificname']
                                row_data['scientificNameID'] = v_match[0]['lsid']
                                row_data.update(dict(zip(w_ranks, [v_match[0].get(key) for key in w_ranks])))
                                row_data.update({'taxonRank': v_match[0]['rank']})
                                #print(old_name+": Unaccepted, using "+valid_name+", "+rank+" ")
                            #else:
                                #print(old_name+": Unaccepted, no valid name, "+rank+" ")
        matches += [row_data]
    matches = pd.DataFrame.from_dict(matches)
    queue.put(matches)


def get_worms_from_aphiaid_or_name_parallel(tax_df, worms_dict, ordered_rank_columns, full_tax_column="taxonomy", like=False, marine_only=False, full_tax_vI=False, n_proc=0):
    
    import multiprocess as mp
    import pandas as pd
    import os
    
    # Use spawn method on Windows to avoid issues with forking
    if os.name == 'nt':  # Windows
        mp.set_start_method('spawn', force=True)
    
    if n_proc == 0:
        # create as many processes as there are CPUs on your machine
        num_processes = mp.cpu_count()
    else:
        num_processes = n_proc
    
    # calculate the chunk size as an integer
    chunk_size = max(1, round(tax_df.shape[0]/num_processes))
    
    # Create chunks
    chunks = [tax_df.iloc[tax_df.index[i:i + chunk_size]] for i in range(0, tax_df.shape[0], chunk_size)]
    
    # Create a manager to handle the queue
    manager = mp.Manager()
    queue = manager.Queue()
    
    # Create and start processes
    procs = []
    for chunk in chunks:
        proc = mp.Process(
            target=get_worms_from_aphiaid_or_name,
            args=(chunk, worms_dict, ordered_rank_columns, queue, full_tax_column, like, marine_only, full_tax_vI)
        )
        procs.append(proc)
        proc.start()
    
    # Collect results from the queue
    results = []
    for _ in range(len(procs)):
        try:
            results.append(queue.get(timeout=600))  # 10-minute timeout
        except Exception as e:
            print(f"Error getting results from queue: {e}")
    
    # Wait for all processes to finish
    for proc in procs:
        proc.join(timeout=10)  # 10-second timeout for joining
        if proc.is_alive():
            print(f"Warning: Process {proc.pid} did not terminate properly, forcing termination")
            proc.terminate()
    
    # Combine results
    new_df = pd.DataFrame()
    for result in results:
        new_df = pd.concat([new_df, result])
    
    return new_df


if __name__ == '__main__':
    pass

