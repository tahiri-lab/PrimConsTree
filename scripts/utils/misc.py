import os


def create_unique_file(base_filename):
    """
    100% ChatGpt here to be honest, don't look very good but do the job
    Take a path, create directories if necessary. If the file already exists,
    create a new unique filename.
    """
    dir_name = os.path.dirname(base_filename)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    # Check if the file already exists
    if not os.path.exists(base_filename):
        # If the file does not exist, create it
        with open(base_filename, "w") as file:
            pass
        return base_filename

    # If the file exists, append an integer to the filename
    file_extension = ""
    if "." in base_filename:
        base_name, file_extension = base_filename.rsplit(".", 1)
        base_filename = base_name
        file_extension = "." + file_extension

    counter = 1
    new_filename = f"{base_filename}_{counter}{file_extension}"

    while os.path.exists(new_filename):
        counter += 1
        new_filename = f"{base_filename}_{counter}{file_extension}"

    # Create the new file with a unique name
    with open(new_filename, "w") as file:
        pass

    return new_filename


def get_alg_id(alg: tuple[str, dict]) -> str:
    """Create an id based on alg additional parameters"""
    return alg[0] + "_" + "_".join([str(k) + str(v) for k, v in alg[1].items()])
