import os
import shutil
import time
import glob


def delete_files_in_subfolder(subfolder, files_to_delete):

    # Delete files in the current subfolder
    for file_to_delete in files_to_delete:
        file_path = os.path.join(subfolder, file_to_delete)
        for filename in glob.glob(file_path):
            try:
                os.remove(filename)
                print(f"Deleted {filename}")
            except FileNotFoundError:
                print(f"{filename} not found")


def get_processor_folders(parent_folder):
    # Get all subfolders with names starting with "processor" followed by an integer
    processor_folders = [
        os.path.join(parent_folder, folder)
        for folder in os.listdir(parent_folder)
        if folder.startswith("processor") and folder[9:].isdigit()
    ]
    return processor_folders


def monitor_processor_folders(parent_folder, files_to_delete, duration_hours,
                              cleaning_interval):

    end_time = time.time() + duration_hours * 3600
    existing_subfolders = {
        folder: set()
        for folder in get_processor_folders(parent_folder)
    }

    while time.time() < end_time and len(existing_subfolders) == 0:

        print('time', end_time - time.time())

        existing_subfolders = {
            folder: set()
            for folder in get_processor_folders(parent_folder)
        }

        time.sleep(cleaning_interval)

    while time.time() < end_time:
        print('time', end_time - time.time())
        # print('existing_subfolders',existing_subfolders)

        for processor_folder in get_processor_folders(parent_folder):
            current_subfolders = set([
                subfolder for subfolder in os.listdir(processor_folder)
                if subfolder.replace('.', '').isnumeric()
            ])

            latest_subfolder = max(current_subfolders, default="")
            current_subfolders.discard(latest_subfolder)

            new_subfolders = current_subfolders - \
                existing_subfolders[processor_folder]

            if len(new_subfolders) > 0:

                for subfolder in new_subfolders:
                    subfolder_path = os.path.join(processor_folder, subfolder)
                    delete_files_in_subfolder(subfolder_path, files_to_delete)

            existing_subfolders[processor_folder] = current_subfolders

        time.sleep(cleaning_interval)


if __name__ == "__main__":
    # Set the path to the parent folder containing "processor" folders
    parent_folder = "./"

    # List of files to be deleted
    files_to_delete = [
        "kineticTheory*", "alphaRhoPhi*", "nut.*", "phi.*", "Theta.*",
        "alphat.*"
    ]

    # Set the duration to run the script in hours
    duration_hours = 0.2

    # Set the cleaning interval in seconds
    cleaning_interval = 5

    # Start monitoring the "processor" folders for new subfolders for the specified duration
    monitor_processor_folders(parent_folder, files_to_delete, duration_hours,
                              cleaning_interval)
