#!/usr/bin/env python
#!/usr/bin/env python
#
# Author   : Antonio M. Ferreira, PhD
#            Center for Simulation and Modeling
#            University of Pittsburgh
# Date     : 22 Aug 14
# Modified : 15 Jan 15
# Now moved to github
#
# Things to do:
#
#  1) In direct mode, make limiting the download speed to disk speed optional
#  2) Make the disk speed check a function and call it periodically to adjust
#     for I/O load
#  3) Build in logic for handling more than one file per CGHub UUID
#

import time
import sys
import os
import getopt
import subprocess
from datetime import datetime

verbose     = 0
gt_debug    = 0
Bytes2MB    = 1048576
Bytes2GB    = 1073741824
TimeFormat  = '%m/%d/%y %I:%M %p'
final_dest  = "/supercell/tcga"
exit_code   = 0
debug       = 1
direct_mode = 1

# Default values
GeneTorrentExecutable = "/usr/bin/gtdownload"
CGQueryExecutable     = "/usr/bin/cgquery"
CredentialFile        = "/arc/users/omedvede/bams_private/cghub.key"
# These values were optimized for firehose6 @ PSC
MAX_WAIT              = 10
num_children          = 24
max_bandwidth         = 10000
local_dir             = "/tmp"
request_file_dir      = "/supercell/bam_requests"
request_file_name     = "bam_status.tsv"
########################################################
# For debugging we'll use local files and directories
#
if debug:
    # For running on firehose6 at PSC
    # CredentialFile    = "/arc/users/aferreir/BAM_download/Olgas.key"
    # final_dest        = "/supercell2/users/aferreir"

    # For running on corelli.sam.pitt.edu
    CredentialFile    = "/home/tony/Documents/Research/PGRR/BAM_Downloads/Olgas.key"
    final_dest        = "/home/tony/Documents/Research/PGRR/BAM_Downloads/DownloadDir"
    num_children      = 6

    request_file_dir  = "."
    request_file_name = "test.tsv"
    verbose           = 1
########################################################

# Parse the command line options
options, trailing_opts = getopt.getopt(sys.argv[1:],
                                      'e:q:C:w:t:b:d:l:n:vhD', [])


class BAMinfo:
    def __init__(self,filename,bamID,filesize,bamsum,start,end,stat,disease_str,barcode_str,library_type,platform_name):
        self.name        = filename
        self.uuid        = bamID
        self.size        = filesize
        self.checksum    = bamsum
        self.start_time  = start
        self.end_time    = end
        self.status      = stat
        self.disease     = disease_str
        self.barcode     = barcode_str
        self.library     = library_type
        self.platform    = platform_name
        self.localname   = disease_str.lower() + '/' + '-'.join(barcode.split('-')[0:3]) + '/' \
                           + '-'.join(barcode.split('-')[0:4])[:-1] + '/' + library_type + '/' \
                           + "CGHub_" + platform_name + '/' + bamID + '/' + filename
    def bamprint(self):
        print ("\n___________________________________________________________________________________")
        print ("BAM info for %s") % self.name
        print ("CGHub UUID        = %s")     % self.uuid
        print ("Size (Bytes)      = %d")     % self.size
        print ("Disease type      = %s")     % self.disease
        print ("Barcode           = %s")     % self.barcode
        print ("Library type      = %s")     % self.library
        print ("Analysis platform = %s")     % self.platform
        print ("Local copy        = .../%s") % self.localname
        print ("Download start    = %s")     % self.start_time
        print ("Download end      = %s")     % self.end_time
        print ("Status            = %s")     % self.status
        print ("___________________________________________________________________________________")
        print ("\n")


def usage():
    print "\nUsage: GT_Download.py [OPTION] [REQUESTS_FILE]"
    print "\nThis script parses a tab-delimited REQUESTS_FILE containing a list of files to be"
    print "downloaded from CGHub, downloads each file to temporary storage, copies it to a"
    print "destination directory specified by REQUESTS_FILE, and then updates the"
    print "REQUESTS_FILE.\n"
    print "  REQUESTS_FILE    is an optional filename containing the download requests (default = %s)" % \
          request_file_name
    print "  -e NAME          full path to alternate gtdownload executable (default = %s)" % GeneTorrentExecutable
    print "  -q NAME          full path to alternate cgquery executable (default = %s)" % CGQueryExecutable
    print "  -C FILE          Credential file to use for CGHub access (default = %s)" % CredentialFile
    print "  -w TIME          wait TIME (in minutes) before assuming a download has stalled (default = %d)" % MAX_WAIT
    print "  -n NUM           Use NUM threads when downloading (default = %d)" % num_children
    print "  -b NUM           Limit download bandwidth to NUM Mb/s (default = %d)" % max_bandwidth
    print "  -d DIR           Look for request file in DIR (default = %s)" % request_file_dir
    print "  -l DIR           local directory to use for initial download (default = %s)" % local_dir
    print "  -t DIR           Target DIR for downloaded files (default = %s)" % final_dest
    print "  -D               Turn off DIRECT download mode (i.e. - cache to local directory first)"
    print "  -v               Verbose mode"
    print "  -h               Show this help message\n\n"



def UpdateRequestsFile(filename,bamname,num_cols,col_index,newdata):
    # A function for updating the .tsv file
    import tempfile,shutil

    BUFFER_SIZE = 10485760 # 10 MB should be enough (for now)

    tmp = tempfile.SpooledTemporaryFile(bufsize=BUFFER_SIZE)
    f = open(filename,'r')
    f.seek(0)
    tmp.seek(0)
    for line in f:
        data = line.strip().split('\t')
        if len(data) < num_cols:
            for i in range(len(data),(num_cols+1)):
                data.append("")
        if data[13] == bamname:
            data[col_index] = newdata
            line = '\t'.join(data)
            line = line + '\n'
        tmp.write(line)
    f.close()
    shutil.move(filename, os.path.join(os.path.dirname(filename),"." + os.path.basename(filename)))
    tmp.seek(0)
    f = open(filename,'w')
    for line in tmp:
        f.write(line)
    f.close()
    # Use flush + fsync to ensure that we have committed the update to disk
    f.flush()
    os.fsync(f.fileno())
    tmp.close()

# Get the start time for the script
script_start_time = datetime.now()

# Process the commandline arguments
if len(trailing_opts) > 0:
    request_file_name = trailing_opts[0]
for opt, arg in options:
    if opt == '-e':
        GeneTorrentExecutable = os.path.abspath(arg)
    elif opt == '-q':
        CGQueryExecutable = os.path.abspath(arg)
    elif opt == '-C':
        CredentialFile = os.path.abspath(arg)
    elif opt == '-w':
        MAX_WAIT = int(arg)
    elif opt == '-n':
        num_children = int(arg)
    elif opt == '-b':
        max_bandwidth = int(arg)
    elif opt == '-d':
        request_file_dir = os.path.abspath(arg)
    elif opt == '-l':
        local_dir = os.path.abspath(arg)
    elif opt == '-t':
        final_dest = os.path.abspath(arg)
    elif opt == '-D':
        direct_mode = 0
    elif opt == '-v':
        verbose = 1
    elif opt == '-h':
        usage()
        quit()
    else:
        print "Unknown option: %s" % opt
        quit()

# Normalize the path to the requests file
RequestsFileName = os.path.abspath(request_file_dir + '/' + request_file_name)

# Dump out the run options for verification
print ("====================================================================================")
print ("=                         GeneTorrent download parameters                          =")
print ("====================================================================================")
print ("GeneTorrent executable      = %s") % GeneTorrentExecutable
print ("CGQuery executable          = %s") % CGQueryExecutable
print ("Credential file             = %s") % CredentialFile
print ("BAM requests file           = %s") % RequestsFileName
if not direct_mode:
    print ("Local directory for caching = %s") % local_dir
print ("Target directory            = %s") % final_dest
print ("Number of child processes   = %d") % num_children
print ("Maximum network bandwidth   = %d MB/s") % max_bandwidth
if direct_mode:
    print ("Data transfer timeout       = %d min") % MAX_WAIT
    print ("Running in direct mode\n")
else:
    print ("Data transfer timeout       = %d min\n") % MAX_WAIT

# Make sure the GeneTorrent executables are where we think they are
if not(os.path.isfile(GeneTorrentExecutable)):
    print ("%s is not a file!") % GeneTorrentExecutable
    exit_code = 10
    sys.exit(exit_code)
if not(os.path.isfile(CGQueryExecutable)):
    print ("%s is not a file!") % CGQueryExecutable
    exit_code = 10
    sys.exit(exit_code)

# Make sure the requests file exists and we can write to it
if not(os.path.isfile(RequestsFileName)):
    print ("%s is not a file!") % RequestsFileName
    exit_code = 10
    sys.exit(exit_code)
elif not(os.access(RequestsFileName, os.W_OK)):
    print ("You do not have permission to write to %s") % RequestsFileName
    exit_code = 10
    sys.exit(exit_code)

# Check that we can write to the cache directory
if not(os.access(local_dir, os.W_OK)) :
    print ("You do not have permission to write to %s") % \
          local_dir
    exit_code = 10
    sys.exit(exit_code)

# Check that we can write to the final destination directory
if not(os.access(final_dest, os.W_OK)):
    print ("You do not have permission to write to %s") % \
          final_dest
    exit_code = 10
    sys.exit(exit_code)

# If we're running in direct mode, test the underlying filesystem and reduce the
# download speed to match the filesystem
if direct_mode:
    if verbose:
        print ("Running filesystem speed test for target directory (%s)") % final_dest
        sys.stdout.flush()
    speedcmd = "dd if=/dev/zero of=%s/GT_Download.speedtest bs=4k count=75000 conv=fdatasync" % final_dest
    if debug:
        print ("DEBUG: speedcmd = %s") % speedcmd
        sys.stdout.flush()
    cmd = subprocess.Popen(speedcmd, shell=True, bufsize=1,\
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    speedout,speederr = cmd.communicate()
    os.remove(final_dest + "/GT_Download.speedtest")
    if cmd.returncode != 0:
        print ("ERROR: Failed on shell command '%s'") % mkdircmd
        sys.stdout.flush()
    else:
        test_bandwidth = int(float(speederr.split()[len(speederr.split())-2]))
        speed_units = speederr.split()[len(speederr.split())-1]
        if speed_units == "GB/s":
            test_bandwidth = test_bandwidth * 1000
        if debug:
            print ("DEBUG: speederr = %s") % speederr
            print ("DEBUG: filesystem bandwidth = %s") % test_bandwidth
        if test_bandwidth < max_bandwidth:
            if verbose:
                print ("Filesystem bandwidth is only %s MB/s. Lowering CGHub download bandwidth to match.") \
                      % test_bandwidth
            max_bandwidth = test_bandwidth


# Parse the requests file and build a list of the downloads to perform
SourceList = list()
RequestsFile = open(RequestsFileName,'r')
RequestsFile.seek(0)
column_names = RequestsFile.readline().strip().split('\t')
num_columns = len(column_names)
for line in RequestsFile:
    data = line.strip().split('\t')
    # Fix data for truncated lines
    if len(data) < num_columns:
        for i in range(len(data),(num_columns+1)):
            data.append("")
    filename = data[column_names.index('filename')]
    BAM_UUID = data[column_names.index('analysis_id')]
    filesize = int(float(data[column_names.index('files_size')]))
    bamsum = data[column_names.index('checksum')]
    stime  = data[column_names.index('start_time')]
    if stime != '':
        stime = datetime.strptime(stime,TimeFormat)
    else:
        stime = datetime.now()
    etime  = data[column_names.index('end_time')]
    if etime != '':
        etime = datetime.strptime(etime,TimeFormat)
    else:
        etime = datetime.now()
    stat     = data[column_names.index('status')]
    disease  = data[column_names.index('disease')]
    barcode  = data[column_names.index('barcode')]
    library  = data[column_names.index('library_type')]
    platform = data[column_names.index('platform_name')]
    current  = BAMinfo(filename,BAM_UUID,filesize,bamsum,stime,etime,stat,disease,barcode,library,platform)
    if current.status == "Finished":
        if verbose:
            print "%-51s: downloaded on %s" % (current.name,str(current.end_time))
    elif current.status == "Live":
        if verbose:
            print "%-51s: downloaded on %s and is now live." % (current.name,str(current.end_time))
    elif current.status == "":
        current.status = "Unknown"
        SourceList.append(current)
    else:
        SourceList.append(current)
RequestsFile.close()

# Dump out the list of downloads to perform
print ("\n\n                                 FILES TO DOWNLOAD")
print ("%-71s %12s") % ('Name','Size (MB)')
print ("====================================================================================")
total_download_size = 0
num_bams = len(SourceList)
for bam in SourceList:
    print ("%-71s %12d") % (bam.name,float(bam.size/Bytes2MB))
    total_download_size += bam.size
print ("\nTotal anticipated download size = %.2f GB (%d UUIDs in all)\n") % \
      (float(total_download_size)/float(Bytes2GB),num_bams)
total_download_size = 0

# Loop through the list of files to download
bam_count = 0
for bam in SourceList:
    bam_count += 1
    source_uuid = bam.uuid
    # Where is this thing finally going to end up?
    final_location = os.path.join(final_dest,bam.localname)

    # Create the gtdownload command
    gt_command = "%s" % GeneTorrentExecutable
    gt_command += " -t"             # Timestamp log messages
    if verbose:
        gt_command += " -v"
    if gt_debug:
        gt_command += " -l stdout:full" # Full logging to stdout
        gt_command += " -vv"            # Detailed progress information
    gt_command += " --max-children %d" % num_children
    gt_command += " --rate-limit %d" % max_bandwidth
    gt_command += " -k %d" % MAX_WAIT
    gt_command += " -c %s" % CredentialFile
    if direct_mode:
        direct_mode_path = os.path.dirname(bam.localname)
        direct_mode_path = os.path.join(os.path.split(direct_mode_path)[:1])[0]
        direct_mode_path = os.path.join(final_dest,direct_mode_path)
        if debug:
            print "DEBUG: direct_mode_path = ",direct_mode_path
        gt_command += " -p %s" % direct_mode_path
        final_location = direct_mode_path + "/" + bam.name
    else:
        gt_command += " -p %s" % local_dir
    gt_command += " -d %s" % source_uuid

    if debug:
        print ("DEBUG: final_location = %s") % final_location

    # Print the run parameters
    print ("\nCurrent download (%d of %d): %s") % (bam_count,num_bams,bam.name)
    print ("====================================================================================")
    print ("CGHub data source = %s") % source_uuid
    if direct_mode:
        print ("Destination       = %s") % (os.path.dirname(final_location) + "/" + bam.uuid)
    else:
        print ("Destination       = %s") % os.path.dirname(final_location)
    sys.stdout.flush()

    # Check to see if the target directory exists and create if necessary
    if not(os.path.exists(os.path.dirname(final_location))):
        if verbose:
            print ("Creating directory %s") % os.path.dirname(final_location)
            sys.stdout.flush()
        mkdircmd = "mkdir -p %s" % os.path.dirname(final_location)
        cmd = subprocess.Popen(mkdircmd, shell=True)
        cmd.wait()
        if cmd.returncode != 0:
            print ("ERROR: Failed on shell command '%s'") % mkdircmd
            sys.stdout.flush()
            # sys.exit(cmd.returncode)
    else:
        if verbose:
            print ("Target directory (%s) exists.") % os.path.dirname(final_location)
            sys.stdout.flush()

    start_time = datetime.now()
    do_download=1
    attempt=0
    MAX_ATTEMPTS=5

    # If the status is Cached, Staged, Finished, or Live, skip this step
    if (bam.status == "Cached" or bam.status == "Staged" or bam.status == "Finished" or bam.status == "Live"):
        print ("This file has already been downloaded. Status = %s") % bam.status
        do_download = 0

    # Check that the source is downloadable on the CGHub side
    cgquery_cmd = "%s \"analysis_id=%s\" -a" % (CGQueryExecutable,bam.uuid)
    if debug:
        print ("\n>>>>>> cgquery command = %s") % cgquery_cmd
        sys.stdout.flush()
    cgquery_process = subprocess.Popen(cgquery_cmd, shell=True, bufsize=1,\
                                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    cgquery_out,cgquery_err = cgquery_process.communicate()
    if gt_debug:
        print (">>>>>> cgquery stdout:\n%s\n\ncgquery stderr:\n%s") % (cgquery_out,cgquery_err)
    if cgquery_process.returncode != 0:
        print ("ERROR: Failed on shell command '%s'") % cgquery_cmd
        sys.stdout.flush()
        # sys.exit(cgquery_process.returncode)
    else:
        # Process the output here
        success_string = "All matching objects are in a downloadable state."
        if (cgquery_out.find(success_string) == -1):
            print ("UUID %s is not in a downloadable state. Skipping this entry.") % bam.uuid
            state_loc = cgquery_out.find("state_count")
            bam.status = cgquery_out[state_loc:].split()[1]
            if debug:
                print ("DEBUG: bam.status = %s") % bam.status
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('status'),bam.status)
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('state'),bam.status)
            do_download = 0
        else:
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('state'),"Live")
            if verbose:
                print ("UUID %s is in a downloadable state.") % bam.uuid

    # Dump some debugging output
    if debug:
        print ("DEBUG: Info after cgquery check")
        bam.bamprint()
        print ("DEBUG: do_download = %d") % do_download
        sys.stdout.flush()

    # Start the actual download here.  The do_download variable is used to
    # short-circuit the download process if necessary (see above)
    while (do_download and attempt < MAX_ATTEMPTS):
        if debug:
            print ("\nDEBUG: gtdownload command = %s\n") % gt_command
            sys.stdout.flush()

        attempt += 1
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('download_attempt_num'),\
                           str(attempt))
        bam.status = "InProcess"
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('status'),bam.status)
        bam.start_time = datetime.now()
        # A little user output for the log file
        print (" --- Starting download at %s") % (bam.start_time.strftime(TimeFormat))
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('start_time'),\
                           bam.start_time.strftime(TimeFormat))
        gt_process = subprocess.Popen(gt_command, shell=True, bufsize=1,\
                                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if verbose:
            print (">>>>>> %s output:") % GeneTorrentExecutable
            # A little user output for the log file
        while (gt_process.poll()  == None):
            out = gt_process.stdout.readline()
            if out == '':
                if gt_process.poll() != None:
                    if direct_mode:
                        cached_name = final_dest
                    else:
                        cached_name = "%s/%s/%s" % (local_dir,bam.uuid,bam.name)
                    if (gt_process.returncode == 0 and os.path.exists(cached_name)):
                        do_download = 0
                        bam.end_time = datetime.now()
                        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                                           column_names.index('end_time'),\
                                           bam.end_time.strftime(TimeFormat))
                        if direct_mode:
                            bam.size = os.path.getsize(os.path.dirname(final_location) + \
                                                       "/" + bam.uuid + "/" + bam.name)
                            total_download_size += bam.size
                            if debug:
                                print ("DEBUG: Getting size of %s (%d)") % (os.path.dirname(final_location) + \
                                                       "/" + bam.uuid + "/" + bam.name, bam.size)
                        else:
                            if debug:
                                print ("DEBUG: Getting size of %s") % cached_name
                            bam.size = os.path.getsize(cached_name)
                            total_download_size += bam.size
                        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                                           column_names.index('files_size'),\
                                           str(bam.size))
                        if direct_mode:
                            bam.status = "Finished"
                        else:
                            bam.status = "Cached"
                        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                                           column_names.index('status'),\
                                           bam.status)
                        break
                    else:
                        print ("\nERROR: Download process failed with exit code %d.  Retrying.  (Attempt %d of %d)\n\n") % \
                              (gt_process.returncode,(attempt+1),MAX_ATTEMPTS)
                        bam.end_time = datetime.now()
                        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                                           column_names.index('end_time'),\
                                           bam.end_time.strftime(TimeFormat))
                        bam.status = "Failed"
                        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                                           column_names.index('status'),\
                                           bam.status)
            else:
                sys.stdout.write(">>>>>> " + out)
                sys.stdout.flush()

    # A little user output for the log file
    print (" --- Finished download at %s") % (bam.end_time.strftime(TimeFormat))

    # How long did this take?
    elapsed = bam.end_time - bam.start_time
    elapsed_time = elapsed.seconds
    (elapsed_days,elapsed_hours) = divmod(elapsed_time,(3600 * 24))
    (elapsed_hours,elapsed_minutes) = divmod(elapsed_hours,3600)
    (elapsed_minutes,elapsed_seconds) = divmod(elapsed_minutes,60)
    print (" --- Elapsed time for download operation = %d days %2d hours %2d minutes %4.1f seconds") % \
          (elapsed_days,elapsed_hours,elapsed_minutes,elapsed_seconds)
    sys.stdout.flush()

    # Clean up the .gto file
    if not (bam.status == "Failed" or bam.status == "suppressed"):
        if direct_mode:
            os.remove(direct_mode_path + "/" + bam.uuid + ".gto")
        else:
            os.remove(local_dir + "/" + bam.uuid + ".gto")

    # What was our data rate?
    if elapsed_time == 0:
        data_rate = 0
    else:
        data_rate = float(bam.size/Bytes2MB) / (bam.end_time - bam.start_time).seconds
    if verbose:
        print (" --- Calculated data rate = %.1f MB/s") % data_rate
    UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                       column_names.index('overall_rate_(MB/s)'),str(data_rate))

    if (bam.status == "Cached" and not direct_mode):
        if verbose:
            print ("Downloaded files(s) will now be copied to %s") % final_location

        # Copy the file and to its final resting place
        copy_start = datetime.now()

        # Copy the files
        # Let's just copy the entire directory instead.  That way we don't have to worry about
        # specific files inside it.  This also makes using rsync more reasonable
        #
        copycmd = "rsync -rt"
        if verbose:
            copycmd += "v"
        copycmd += " %s/%s %s" % (local_dir,bam.uuid,\
                                 '/'.join(os.path.split(os.path.dirname(final_location))[:-1]))
        if debug:
            print ("DEBUG: Copy command = %s") % copycmd
            sys.stdout.flush()
        cmd = subprocess.Popen(copycmd, shell=True)
        cmd.wait()
        copy_end = datetime.now()
        if cmd.returncode != 0:
            print ("ERROR: Failed on shell command '%s'") % copycmd
            sys.stdout.flush()
        bam.status = "Staged"
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('status'),bam.status)
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('pgrr_file_path'),\
                           os.path.dirname(bam.localname))
        elapsed_copy = copy_end - copy_start
        elapsed_time = elapsed_copy.seconds
        (elapsed_days,elapsed_hours) = divmod(elapsed_time,(3600 * 24))
        (elapsed_hours,elapsed_minutes) = divmod(elapsed_hours,3600)
        (elapsed_minutes,elapsed_seconds) = divmod(elapsed_minutes,60)
        print (" --- Elapsed time for copy operation     = %d days %2d hours %2d minutes %4.1f seconds") % \
              (elapsed_days,elapsed_hours,elapsed_minutes,elapsed_seconds)
        sys.stdout.flush()
    elif bam.status == "suppressed" and verbose:
        print ("This file is not currenlty available from CGHub. Status = %s") % bam.status
    elif (bam.status == "Finished" or bam.status == "Staged") and verbose:
        print ("This file is already in the final destination directory. Status = %s") % bam.status
    elif verbose:
        print ("This location/availability of this file is unknown. Status = %s") % bam.status

    if direct_mode:
        UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                           column_names.index('pgrr_file_path'),\
                           os.path.dirname(bam.localname))
    
    if (bam.status == "Staged" and not direct_mode):
        # Perform an md5sum on the file and compare to BAMinfo
        if verbose:
            print ("Checking md5sum of file after copy.")
            sys.stdout.flush()
        md5_start = datetime.now()
        md5cmd = "md5sum %s" % final_location
        md5_process = subprocess.Popen(md5cmd, shell=True, bufsize=1,\
                                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = md5_process.communicate()
        if md5_process.returncode != 0:
            print ("ERROR: Failed on shell command '%s'") % md5cmd
            my_md5 = 0
            sys.stdout.flush()
            # sys.exit(md5_process.returncode)
        else:
            my_md5 = out.split()[0]
        md5_end = datetime.now()
        elapsed_md5 = md5_end - md5_start
        elapsed_time = elapsed_md5.seconds
        (elapsed_days,elapsed_hours) = divmod(elapsed_time,(3600 * 24))
        (elapsed_hours,elapsed_minutes) = divmod(elapsed_hours,3600)
        (elapsed_minutes,elapsed_seconds) = divmod(elapsed_minutes,60)
        print (" --- Elapsed time for MD5 checksum       = %d days %2d hours %2d minutes %4.1f seconds") % \
              (elapsed_days,elapsed_hours,elapsed_minutes,elapsed_seconds)
        sys.stdout.flush()

        if my_md5 != bam.checksum:
            print ("ERROR: Checksum does not match!")
            print ("ERROR:   - Reference  = %s") % bam.checksum
            print ("ERROR:   - Calculated = %s") % my_md5
            bam.status = "Failed"
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('status'),bam.status)
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('end_time'),"")
            exit_code = 5
            sys.stdout.flush()
            # sys.exit(exit_code)
        else:
            bam.status = "Finished"
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('status'),bam.status)
            data_rate = float((bam.size / (md5_end - start_time).seconds)/Bytes2MB)
            download_speed = "%.2f" % data_rate
            UpdateRequestsFile(RequestsFileName,bam.name,num_columns,\
                               column_names.index('overall_rate_(MB/s)'),\
                               download_speed)
            if verbose:
                print ("Checksum passed.  Removing cached file.")
                sys.stdout.flush()
            # Remove the cached copy of the file
            tempcmd = "rm -rf %s/%s*" % (local_dir,source_uuid)
            temp = subprocess.Popen(tempcmd, shell=True)
            temp.wait()
            if temp.returncode != 0:
                print ("ERROR: Failed on shell command '%s'") % tempcmd

            if verbose:
                bam.bamprint()

        end_time = datetime.now()
        elapsed_time = (end_time - start_time).seconds
        (elapsed_days,elapsed_hours) = divmod(elapsed_time,(3600 * 24))
        (elapsed_hours,elapsed_minutes) = divmod(elapsed_hours,3600)
        (elapsed_minutes,elapsed_seconds) = divmod(elapsed_minutes,60)
        print (" --- Total time for dataset              = %d days %2d hours %2d minutes %4.1f seconds") % \
              (elapsed_days,elapsed_hours,elapsed_minutes,elapsed_seconds)
        if bam.status != "Finished":
            effective_data_rate = 0
        else:
            effective_data_rate = float(bam.size/Bytes2MB) / elapsed_time
        print (" --- Effective data rate = %.1f MB/s") % data_rate
        print ("====================================================================================\n\n")
        sys.stdout.flush()

    # Prune the directory tree to remove any empty directories if something failed
    # NOTE: Olga requested this behavior
    if not (bam.status == "Finished" or bam.status == "Staged"):
        if debug:
            print ("DEBUG: Problem with download of %s (status=%s)") % (bam.uuid,bam.status)
            print ("DEBUG: Checking for empty directories in %s") % final_dest
        done = 0
        while not done:
            empty_dirs = list()
            dirinfo = os.walk(final_dest)
            for current_dir in dirinfo:
                numdirs  = len(current_dir[1])
                numfiles = len(current_dir[2])
                if numdirs == 0 and numfiles == 0:
                    empty_dirs.append(current_dir[0])
            if len(empty_dirs) == 0:
                done = 1
            else:
                if debug:
                    print ("DEBUG: Empty directory list:")
                    for i in empty_dirs:
                        print ("DEBUG: ----> %s") % i
                for empty_path in empty_dirs:
                    import shutil
                    if debug:
                        print ("DEBUG: Removing empty directory %s") % empty_path
                    shutil.rmtree(empty_path)

# Print some summary output
StatusList = list()
for bam in SourceList:
    included = 0
    for i in xrange(len(StatusList)):
        name,count = StatusList[i]
        if bam.status == name:
            StatusList[i][1] = count + 1
            included = 1
    if not included:
        StatusList.append([bam.status,1])

print ("\n\n")
print ("    DOWNLOAD SUMMARY    ")
print ("BAM status         Count")
print ("------------------------")
for name,count in StatusList:
    print ("%-15s  %7d") % (name,count)
print ("------------------------")
print ("%.2f GB total") % (float(total_download_size)/float(Bytes2GB))

end_time = datetime.now()
elapsed_time = (end_time - script_start_time).seconds
(elapsed_days,elapsed_hours) = divmod(elapsed_time,(3600 * 24))
(elapsed_hours,elapsed_minutes) = divmod(elapsed_hours,3600)
(elapsed_minutes,elapsed_seconds) = divmod(elapsed_minutes,60)
print ("\nTotal walltime for GT_Download           = %d days %2d hours %2d minutes %4.1f seconds") % \
      (elapsed_days,elapsed_hours,elapsed_minutes,elapsed_seconds)
print ("====================================================================================\n")
