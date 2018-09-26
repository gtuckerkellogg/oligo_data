import os
import shutil
import tempfile

# import inspect
# import argparse


#################
##   CLASSES   ##
#################

class TempDir:

    def __init__(self, *, path, delete):

        self.path = path
        self.delete = delete

    def __repr__(self):

        return self.path

    def clean_up(self):

        if self.delete:
            shutil.rmtree(self.path)


###################
##   FUNCTIONS   ##
###################

def not_none(*values):

    """
    Accepts arbitrary number of values. Returns index of first encountered value that
    is NOT None.
    """
    output = []
    for i, v in enumerate(values):
        if v is not None:
            output.append(i)
    return output


def arg_not_none(func, *values):
    
    """
    Checks for valid mutually exclusive arguments (i.e. that mutually exclusive arguments
    are filled in and that EXACTLY 1 is not None)
    """
    valid_i = not_none(*values)
    if len(valid_i) == 1:
        return valid_i[0]
    
    elif not valid_i:
        raise Exception("No valid input for mutually exclusive arguments obtained for function {}.".format(func.__name__))
    else:
        raise Exception("Unexpected excess of mutually exclusive arguments obtained for function {}.".format(func.__name__))


def join_ele(separator, *ele):

    """
    Accepts any element that produces a valid string output when passed into
    function str() EXCEPT the empty string ''.
    """

    ele = [str(e) for e in ele if str(e) != '']
    return separator.join(ele)


def contains_any(ref, *searches):

    """
    Accepts a reference iterable and any number of searches.
    Returns True if at least one of the searches can be found in the iterable.
    (Not deep search. Only searches first level. E.g. char in str, ele in outermost list)
    """
    for search in searches:
        if search in ref:
            return True
    return False


def invert_dict(dictionary):
    return {} if not dictionary else \
        dict((v, k) for k, v in dictionary.items())


def extract_cols(*cols, fname=None, delimiter='', index_start=0):

    f = open(fname, 'r')
    cols = sorted(cols)
    output = []

    for line in f:
        split_cols = line[:-1].split(delimiter)
        curr_line = []
        for col in cols:
            try:
                curr_line.append(split_cols[col-index_start])
            except KeyError:
                continue
        output.append(curr_line)

    return output


# assumes uniform number of columns per row
def transpose(iterable):

    if not iterable:
        return []

    output = []
    cols = len(iterable[0])

    for col in range(cols):
        output.append([row[col] for row in iterable])

    return output


def make_dir(*, dir_path, path_exists_message='', dir_desc="directory", check=True):

    """
    Generic function to check presence of directory (and confirm user choice if it exists)
    before allowing it to be made/used. Returns path to directory.

    'path_exists_message' and 'dir_desc' are only used if 'check' is set to True.
    """

    if not dir_path:
        dir_path = os.getcwd()
    
    # request new directory path if check = True and specified path already exists.
    #   Stop only when user explicitly instructs to use the specified directory
    #   OR user chooses a non-existent path
    if check:
        while os.path.exists(os.path.abspath(dir_path)):
            query = input(join_ele('\n', "Directory '{}' already exists.".format(os.path.abspath(dir_path)),
                          path_exists_message,
                          "Press C to continue anyway. Press any other key to specify another directory. "))
            # exit if user explicitly instructs to use existing directory
            if query == 'C' or query == 'c':
                break
            # get new path to directory
            while True:
                dir_path = input("Please provide another path for a {}:\n".format(dir_desc))
                query = input("Confirm {} '{}'? (Y/N) ".format(dir_desc, os.path.abspath(dir_path)))
                if query == 'Y' or query == 'y':
                    break

    # create directory
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path


def make_temp_dir(*, temp_dir, del_dir=True, project_name=''):

    """
    Creates a temporary directory for storing intermediate files and outputs a TempDir
    object. Please include name of temporary directory in 'temp_dir' variable.

    If del_dir is set to 'True', directory will be deleted during clean up.
    """

    # if directory is to be deleted, make temporary directory according to python built-in
    if del_dir:
        return TempDir(path=tempfile.mkdtemp(), delete=True)

    # else manually create directory + check for path validity
    else:

        # if path to temporary directory is not specified OR
        #   if specified path is the working directory, create "temporary" directory
        #   in working directory with default name
        if not temp_dir or (os.path.abspath(os.getcwd()) == os.path.abspath(temp_dir)):
            temp_dir = os.path.join(os.getcwd(),
                                    join_ele('_', project_name, 'temp'))

        # create "temporary" directory
        temp_dir = make_dir(dir_path=temp_dir,
                            path_exists_message="Based on current settings, the temporary directory WILL {}be deleted.".format('' if del_dir else "NOT "),
                            dir_desc="temporary directory")

        # create TempDir object and return it
        return TempDir(path=temp_dir, delete=del_dir)


# removes specific columns from file, writes to another file
def remove_file_col(*, delimiter, fname_in, fname_out, cols_to_remove=(), index_start=0):

    """
    Accepts path to file (fname_in) delimited consistently by 'delimiter' and a path to
    write output file (fname_out), along with a sequence of columns to be removed
    (cols_to_remove). Default behaviour assumes indexing starts from 0 unless specified
    otherwise using keyword argument 'index_start'.
    """

    # open relevant files
    file_in = open(fname_in, 'r')
    file_out = open(fname_out, "w+")

    # iterate through lines in file_in to remove specified columns
    for entry in file_in.readlines():
        cols = remove_indices(entry.split(delimiter), cols_to_remove,
                              index_start=index_start)
        file_out.write(delimiter.join(cols))

    # close files
    file_in.close()
    file_out.close()


# removes elements of list at arbitrarily specified indices
def remove_indices(iterable, indices_to_remove, index_start=0):

    """
    Accepts an iterable that can be indexed by position and a sequence of indices to
    remove, and returns the iterable (order-conserved) excluding elements at specified
    index positions.

    NOTE: The iterable will be converted to a LIST if it is not a string or tuple.
    """

    # check types, convert to list when iterable is not a string
    if not isinstance(iterable, str) and not isinstance(iterable, tuple):
        iterable = list(iterable)

    # remove elements at specified indices
    for i in sorted(map(lambda x: x-index_start, indices_to_remove), reverse=True):
        iterable = iterable[:i] + iterable[i+1:]
    return iterable


def add_at_index(iterable, index_to_add, new_ele, index_start=0):

    """
    Accepts an iterable that can be indexed by position and an index (int) to add
    variable 'new_ele' to, and returns the iterable (order-conserved) with the
    'new_ele' at the specified index. Elements with indices >= 'index_to_add' will
    be shifted the corresponding number of indices (i.e. len(new_ele).

    NOTE: The iterable will be converted to a LIST if it is not a string.

    NOTE: variable 'new_ele' should be of the same type as the iterable. Or, at least
    it should be able to be converted to the iterable type.
    """

    # check types, convert to list where necessary
    try:
        type(iterable)(new_ele)
    except TypeError:
        raise TypeError("Variables 'iterable' ({}) and 'new_ele' ({}) are incompatible.".format(type(iterable), type(new_ele)))
    if not isinstance(iterable, str) and not isinstance(iterable, tuple):
        iterable = list(iterable)

    iterable = iterable[:index_to_add-index_start] + \
               type(iterable)(new_ele) + \
               iterable[index_to_add+index_start+1:]
    return iterable


def replace_at_index(iterable, indices_to_replace, new_ele, index_start=0):

    """
    Accepts an iterable that can be indexed by position and a sequence of indices to
    replace, and returns the iterable (order-conserved) with the elements at
    specified indices replaced.
    """

    iterable = remove_indices(iterable, indices_to_replace, index_start=index_start)
    for i in sorted(map(lambda x: x-index_start, indices_to_replace)):
        iterable = add_at_index(iterable, i, new_ele)
    return iterable


def replace_ele(iterable, ele, new_ele, index_start=0):

    """
    Accepts an iterable that can be indexed by position and an element (variable 'ele')
    to replace, and returns the iterable (order-conserved) with the specified
    elements replaced with 'new_ele'.
    """

    indices_to_replace = [i for i, v in enumerate(iterable) if v == ele]
    iterable = replace_at_index(iterable, indices_to_replace, new_ele, index_start=index_start)
    return iterable
