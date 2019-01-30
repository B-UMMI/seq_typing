import pickle
import traceback
import shlex
import subprocess
from threading import Timer
import shutil
import time
import functools
import os.path
import sys
import argparse


def start_logger(workdir):
    time_str = time.strftime("%Y%m%d-%H%M%S")
    sys.stdout = Logger(workdir, time_str)
    logfile = sys.stdout.getLogFile()
    return logfile, time_str


class Logger(object):
    def __init__(self, out_directory, time_str):
        self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
        self.terminal = sys.stdout
        self.log = open(self.logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        pass

    def getLogFile(self):
        return self.logfile


def checkPrograms(programs_version_dictionary):
    print('\n' + 'Checking dependencies...')
    programs = programs_version_dictionary
    which_program = ['which', '']
    listMissings = []
    for program in programs:
        which_program[1] = program
        run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program, False, None, False)
        if not run_successfully:
            listMissings.append(program + ' not found in PATH.')
        else:
            print(stdout.splitlines()[0])
            if programs[program][0] is None:
                print(program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0])
            else:
                if program.endswith('.jar'):
                    check_version = ['java', '-jar', stdout.splitlines()[0], programs[program][0]]
                    programs[program].append(stdout.splitlines()[0])
                else:
                    check_version = [stdout.splitlines()[0], programs[program][0]]
                run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version, False, None, False)
                if stdout == '':
                    stdout = stderr
                if program == 'wget':
                    version_line = stdout.splitlines()[0].split(' ', 3)[2]
                else:
                    version_line = stdout.splitlines()[0].split(' ')[-1]
                replace_characters = ['"', 'v', 'V', '+']
                for i in replace_characters:
                    version_line = version_line.replace(i, '')
                print(program + ' (' + version_line + ') found')
                if programs[program][1] == '>=':
                    program_found_version = version_line.split('.')
                    program_version_required = programs[program][2].split('.')
                    if len(program_version_required) == 3:
                        if len(program_found_version) == 2:
                            program_found_version.append(0)
                        else:
                            program_found_version[2] = program_found_version[2].split('_')[0]
                    for i in range(0, len(program_version_required)):
                        if int(program_found_version[i]) > int(program_version_required[i]):
                            break
                        elif int(program_found_version[i]) == int(program_version_required[i]):
                            continue
                        else:
                            listMissings.append('It is required ' + program + ' with version ' +
                                                programs[program][1] + ' ' + programs[program][2])
                else:
                    if version_line != programs[program][2]:
                        listMissings.append('It is required ' + program + ' with version ' + programs[program][1] +
                                            ' ' + programs[program][2])
    return listMissings


def required_programs(programs_version_dictionary):
    missing_programs = checkPrograms(programs_version_dictionary)
    if len(missing_programs) > 0:
        sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missing_programs))


def general_information(script_name, logfile, version, outdir, time_str):
    # Check if output directory exists

    print('\n' + '==========> seq_typing <==========')
    print('\n' + 'Program start: ' + time.ctime())

    # Tells where the logfile will be stored
    print('\n' + 'LOGFILE:')
    print(logfile)

    # Print command
    print('\n' + 'COMMAND:')
    script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), script_name)
    print(sys.executable + ' ' + ' '.join(sys.argv))

    # Print directory where programme was lunch
    print('\n' + 'PRESENT DIRECTORY:')
    present_directory = os.path.abspath(os.getcwd())
    print(present_directory)

    # Print program version
    print('\n' + 'VERSION:')
    script_version_git(version=version, current_directory=present_directory, script_path=script_path)

    return script_path


def setPATHvariable(doNotUseProvidedSoftware, script_path):
    path_variable = os.environ['PATH']
    script_folder = os.path.dirname(script_path)
    # Set path to use provided softwares
    if not doNotUseProvidedSoftware:
        bowtie2 = os.path.join(script_folder, 'src', 'bowtie2-2.2.9')
        samtools = os.path.join(script_folder, 'src', 'samtools-1.3.1', 'bin')
        bcftools = os.path.join(script_folder, 'src', 'bcftools-1.3.1', 'bin')

        os.environ['PATH'] = str(':'.join([bowtie2, samtools, bcftools, path_variable]))

    # Print PATH variable
    print('\n' + 'PATH variable:')
    print(os.environ['PATH'])


def script_version_git(version, current_directory, script_path, no_git_info=False):
    """
    Print script version and get GitHub commit information

    Parameters
    ----------
    version : str
        Version of the script, e.g. "4.0"
    current_directory : str
        Path to the directory where the script was start to run
    script_path : str
        Path to the script running
    no_git_info : bool, default False
        True if it is not necessary to retreive the GitHub commit information

    Returns
    -------

    """
    print('Version {}'.format(version))

    if not no_git_info:
        try:
            os.chdir(os.path.dirname(os.path.dirname(script_path)))
            command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
            run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, 15, False)
            print(stdout)
            command = ['git', 'remote', 'show', 'origin']
            run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, 15, False)
            print(stdout)
        except:
            print('HARMLESS WARNING: git command possibly not found. The GitHub repository information will not be'
                  ' obtained.')
        finally:
            os.chdir(current_directory)


def runTime(start_time):
    end_time = time.time()
    time_taken = end_time - start_time
    hours, rest = divmod(time_taken, 3600)
    minutes, seconds = divmod(rest, 60)
    print('Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's')
    return round(time_taken, 2)


def timer(function, name):
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        print('\n' + 'RUNNING {0}\n'.format(name))
        start_time = time.time()

        results = list(function(*args, **kwargs))  # guarantees return is a list to allow .insert()

        time_taken = runTime(start_time)
        print('END {0}'.format(name))

        results.insert(0, time_taken)
        return results
    return wrapper


def removeDirectory(directory):
    if os.path.isdir(directory):
        shutil.rmtree(directory)


def saveVariableToPickle(variableToStore, pickleFile):
    with open(pickleFile, 'wb') as writer:
        pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
    with open(pickleFile, 'rb') as reader:
        variable = pickle.load(reader)
    return variable


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            print('Exception in ' + func.__name__)
            print(e)

            exc_type, exc_value, exc_tb = sys.exc_info()
            print(''.join(traceback.format_exception(exc_type, exc_value, exc_tb)))
    return wrapped_func


def kill_subprocess_Popen(subprocess_Popen, command):
    print('Command run out of time: ' + str(command))
    subprocess_Popen.kill()


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None, print_comand_True):
    run_successfully = False
    if not isinstance(command, str):
        command = ' '.join(command)
    command = shlex.split(command)

    if print_comand_True:
        print('Running: ' + ' '.join(command))

    if shell_True:
        command = ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    not_killed_by_timer = True
    if timeout_sec_None is None:
        stdout, stderr = proc.communicate()
    else:
        time_counter = Timer(timeout_sec_None, kill_subprocess_Popen, args=(proc, command,))
        time_counter.start()
        stdout, stderr = proc.communicate()
        time_counter.cancel()
        not_killed_by_timer = time_counter.isAlive()

    if proc.returncode == 0:
        run_successfully = True
    else:
        if not print_comand_True and not_killed_by_timer:
            print('Running: ' + str(command))
        if len(stdout) > 0:
            print('STDOUT')
            print(stdout.decode("utf-8"))
        if len(stderr) > 0:
            print('STDERR')
            print(stderr.decode("utf-8"))
    return run_successfully, stdout.decode("utf-8"), stderr.decode("utf-8")


def required_length(tuple_length_options, argument_name):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) not in tuple_length_options:
                msg = 'argument {argument_name} requires one of the following number of arguments:' \
                      ' {tuple_length_options}'.format(argument_name=argument_name,
                                                       tuple_length_options=tuple_length_options)
                parser.error(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def arguments_choices_words(argument_choices, argument_name):
    """
    Check if the words passed using the argument, are between those allowed

    Parameters
    ----------
    argument_choices : list
        List with allowed words, e.g. ['escherichia coli', 'streptococcus agalactiae']
    argument_name : str
        Argument name, e.g. '--species'

    Returns
    -------
    ArgumentsChoicesWords : list
        List with strings passed using the argument, e.g. ('escherichia', 'coli')
    """
    class ArgumentsChoicesWords(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            string_value = ' '.join(values)
            if string_value not in argument_choices:
                msg = 'argument {argument_name}: invalid choice: {string_value} (choose from' \
                      ' {argument_choices})'.format(argument_name=argument_name, string_value=string_value,
                                                    argument_choices=argument_choices)
                parser.error(msg)
            setattr(args, self.dest, values)
    return ArgumentsChoicesWords


def get_sequence_information(fasta_file, length_extra_seq):
    sequence_dict = {}
    headers = {}

    with open(fasta_file, 'rtU') as reader:
        blank_line_found = False
        sequence_counter = 0
        temp_sequence_dict = {}
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if not blank_line_found:
                    if line.startswith('>'):
                        if len(temp_sequence_dict) > 0:
                            if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
                                sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
                                headers[temp_sequence_dict.values()[0]['header'].lower()] = sequence_counter
                            else:
                                print(temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0')
                            temp_sequence_dict = {}

                        if line[1:].lower() in headers:
                            sys.exit('Found duplicated sequence headers')

                        sequence_counter += 1
                        temp_sequence_dict[sequence_counter] = {'header': line[1:].lower(), 'sequence': '', 'length': 0}
                    else:
                        temp_sequence_dict[sequence_counter]['sequence'] += line.upper()
                        temp_sequence_dict[sequence_counter]['length'] += len(line)
                else:
                    sys.exit('It was found a blank line between the fasta file above line ' + line)
            else:
                blank_line_found = True

        if len(temp_sequence_dict) > 0:
            if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
                sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
                headers[temp_sequence_dict.values()[0]['header'].lower()] = sequence_counter
            else:
                print(temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0')

    return sequence_dict, headers


def simplify_sequence_dict(sequence_dict):
    simple_sequence_dict = {}
    for counter, info in sequence_dict.items():
        simple_sequence_dict[info['header']] = info
        del simple_sequence_dict[info['header']]['header']
    return simple_sequence_dict


def chunkstring(string, length):
    return (string[0 + i:length + i] for i in range(0, len(string), length))


def clean_headers_sequences(sequence_dict):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    # print 'Checking if reference sequences contain ' + str(problematic_characters) + '\n'

    headers_changed = False
    new_headers = {}
    for i in sequence_dict:
        if any(x in sequence_dict[i]['header'] for x in problematic_characters):
            for x in problematic_characters:
                sequence_dict[i]['header'] = sequence_dict[i]['header'].replace(x, '_')
            headers_changed = True
        new_headers[sequence_dict[i]['header'].lower()] = i

    if headers_changed:
        print('At least one of the those characters was found. Replacing those with _' + '\n')

    return sequence_dict, new_headers


class Bcolors_print(object):
    # The entire table of ANSI color codes - https://gist.github.com/chrisopedia/8754917
    colours = {'BOLD': '\033[1m', 'HEADER': '\033[1;96m', 'OKBLUE': '\033[94m', 'OKGREEN': '\033[92m',
               'WARNING': '\033[93m', 'FAIL': '\033[91m', 'UNDERLINE': '\033[4m', 'ENDC': '\033[0m'}

    def __init__(self, string, color_code):
        # self.string_bcolors = self.colours[color_code] + string + self.colours['ENDC']
        print(self.colours[color_code] + string + self.colours['ENDC'])


def clean_header(header, problematic_characters):
    new_header = header
    if any(x in header for x in problematic_characters):
        for x in problematic_characters:
            new_header = new_header.replace(x, '_')
    return header, new_header


def parse_reference(reference, problematic_characters):
    reference_dict = {}
    headers_correspondence = {}
    with open(reference, 'rtU') as reader:
        header = None
        sequence = ''
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if line.startswith('>'):
                    if header is not None:
                        reference_dict[header] = sequence
                    original_header, new_header = clean_header(line[1:], problematic_characters)
                    if new_header in headers_correspondence:
                        sys.exit('Possible conflicting sequence header in {reference} file:\n'
                                 '{original_header} header might be the same as {first_header} header after problematic'
                                 ' characters ({problematic_characters}) replacement (new header: {new_header})'
                                 ''.format(reference=reference, original_header=original_header,
                                           first_header=headers_correspondence[new_header],
                                           problematic_characters=problematic_characters, new_header=new_header))
                    header = str(new_header)
                    headers_correspondence[header] = str(original_header)
                    sequence = ''
                else:
                    sequence += line.replace(' ', '').upper()
        if len(sequence) > 0:
            reference_dict[header] = sequence
    return reference_dict, headers_correspondence


def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    return ''.join(list(map(lambda base: complement_dict[base], list(seq))))


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    return complement(reverse(seq))
