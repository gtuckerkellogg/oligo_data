import os
import shutil
import tempfile



#################
##   CLASSES   ##
#################

class TempDir():

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

def join_ele(separator, *ele):

    """
    Accepts any element that produces a valid string output when passed into
    function str() EXCEPT the empty string ''.
    """

    ele = [str(e) for e in ele if str(e) != '']
    return separator.join(ele)


def make_dir(*, dir_path, path_exists_message = '', dir_desc = "directory", check = True):

    """
    Generic function to check presence of directory (and confirm user choice if it exists)
    before allowing it to be made/used. Returns path to directory.

    'path_exists_message' and 'dir_desc' are only used if 'check' is set to True.
    """
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


def make_temp_dir(*, temp_dir, del_dir = True, project_name = ''):

    """
    Creates a temporary directory for storing intermediate files and outputs a TempDir
    object. Please include name of temporary directory in 'temp_dir' variable.

    If del_dir is set to 'True', directory will be deleted during clean up.
    """

    # if directory is to be deleted, make temporary directory according to python built-in
    if del_dir:
        return TempDir(tempfile.mkdtemp(), delete = True)

    # else manually create directory + check for path validity
    else:

        # if path to temporary directory is not specified OR
        #   if specified path is the working directory, create "temporary" directory
        #   in working directory with default name
        if not temp_dir or (os.path.abspath(os.getcwd()) == os.path.abspath(temp_dir)):
            temp_dir = os.path.join(os.getcwd(),
                                    join_str('_', project_name, 'temp'))

        # create "temporary" directory
        temp_dir = make_dir(dir_path = temp_dir,
                            path_exists_message = "Based on current settings, the temporary directory WILL {}be deleted.".format('' if del_dir else "NOT "),
                            dir_desc = "temporary directory")

        # create TempDir object and return it
        return TempDir(path = temp_dir, delete = del_dir)


# removes specific columns from file, writes to another file
def remove_file_col(*, delimiter, fname_in, fname_out, cols_to_remove = (), index_start = 0):

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
                              index_start = index_start)
        file_out.write(delimiter.join(cols))

    # close files
    file_in.close()
    file_out.close()


# removes elements of list at arbitrarily specified indices
def remove_indices(iterable, indices_to_remove, index_start = 0):

    """
    Accepts an iterable that can be indexed by position and a sequence of indices to
    remove, and returns a list of elements in the iterable (order-conserved) excluding
    elements at specified index positions.
    """

    iterable = list(iterable)
    for i in sorted(map(lambda x: x-index_start, indices_to_remove), reverse = True):
        del iterable[i]
    return iterable
