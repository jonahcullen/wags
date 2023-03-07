
import re

# https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
def natural_sort(list, key=lambda s:s):
    """
    Sort the list into natural alphanumeric order.
    """
    def get_alphanum_key_func(key):
        convert = lambda text: int(text) if text.isdigit() else text 
        return lambda s: [convert(c) for c in re.split('([0-9]+)', key(s))]
    sort_key = get_alphanum_key_func(key)
    list.sort(key=sort_key)

def collapse_intervals(l, target_length=10e6):
    """
    Collapses list of intervals into sizes less than the target length.
    """
    tmp = []
    full = []
    curr_length = 0

    # chr with single interval - return interval without length
    if len(l) == 1:
        return [l[0][:-1]]

    for i in range(len(l)-1):
        # add first interval start if empty
        if not tmp:
            tmp.append(l[i][0])

        curr_length += l[i][2]
        
        if curr_length < target_length:
            tmp.append(l[i][1])
            
            if l[i+1][2] > target_length:
                # close the current interval
                tmp.append(l[i][1])
                full.append([tmp[0], tmp[-1]])
                # reset
                curr_length = 0
                tmp = []
            
            # if second to last interval, add last and end
            if i+2 == len(l):
                if tmp:
                    tmp.append(l[i+1][1])
                    full.append([tmp[0], tmp[-1]])
                # add last interval
                else:   
                    full.append([l[i+1][0], l[i+1][1]])
                
        elif curr_length > target_length:
            tmp.append(l[i][1])
            full.append([tmp[0], tmp[-1]])
            # reset
            curr_length = 0
            tmp = []
            
            
            if i+2 == len(l):
                # add last interval
                full.append([l[i+1][0], l[i+1][1]])
                
    return full
