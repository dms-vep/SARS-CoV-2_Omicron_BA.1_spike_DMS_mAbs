import datetime
import ftplib
import os
import tarfile
import pandas as pd

def make_tar_file(fastqs, filename):
    """
    Make the compressed *.tar file with runs to upload to SRA.

    Parameters
    ----------
    fastqs: pd.DataFrame
        Table of the fastq file names and paths for upload to the SRA.
    """
    
    # Add all of the fastq files to upload to the tar directory
    try:
        with tarfile.open(filename, mode='w') as f:
            for i, tup in enumerate(fastqs.itertuples()):
                if (i+1) % 10 == 0 or i == 0 or i+1 == len(fastqs):
                    print(f"Adding file {i + 1} of {len(fastqs)} to {filename}")
                f.add(tup.filename_fullpath, arcname=tup.filename)
            print(f"Added all files to {filename}\n")
    except:
        if os.path.isfile(filename):
           raise 
    print(f"The size of {filename} is {os.path.getsize(filename) / 1e9:.1f} GB\n")

    # Check that all of the files are in the tarred directory
    with tarfile.open(filename) as f:
        files_in_tar = set(f.getnames())
    if files_in_tar == set(fastqs['filename']):
        print(f"{filename} contains all {len(files_in_tar)} expected files.\n")
    else:
        raise Exception(f"{filename} does not have all the expected files.")
        
    print("Finished preparing the tar file for upload.")


def upload_via_ftp(tar_path,
                   ftp_username,
                   ftp_account_folder,
                   ftp_subfolder,
                   ftp_address='ftp-private.ncbi.nlm.nih.gov',
                   ftp_password='ftp_password.txt'
):
    """
    Upload the *.tar file of fastqs to the SRA with FTP.

    Parameters
    ----------
    tar_path: str
        Path to the tarred directory of fastqs.
    ftp_username: str
        Username provided by the SRA.
    ftp_account_folder: str
        Account folder provided by the SRA. 
    ftp_subfolder: str
        Unique sub-folder name chosen for upload.
    """
    
    # Get the password 
    with open(ftp_password) as f:
        password = f.read().strip()
        
    # Start the upload - this can take awhile
    print(f"Starting upload at {datetime.datetime.now()}")

    with ftplib.FTP(ftp_address) as ftp:
        ftp.login(user=ftp_username,
                  passwd=password,
                  )
        ftp.cwd(ftp_account_folder)
        ftp.mkd(ftp_subfolder)
        ftp.cwd(ftp_subfolder)
        with open(tar_path, 'rb') as f:
            ftp.storbinary(f"STOR {tar_path}", f)

    print(f"Finished upload at {datetime.datetime.now()}")
