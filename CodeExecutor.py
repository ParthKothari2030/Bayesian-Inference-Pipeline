import os
import subprocess

"""
Module name: Code executor.
Function: Executing the compiled C++/C code via shell process.
Uses: LIM simulator and power spectrum exection for MCMC.
"""

class CodeExecutor:
    """
    CodeExecutor class executes the C++ and C-based compiled code via subprocess (shell) by running the executable file generated
    via "make" library.
    """
    def __init__(self, directory, executable, *params):
        """
        Initialize the CodeExecutor with the directory, executable, and parameters.

        Args:
        directory (str): The directory where the executable is located.
        executable (str): The name of the executable.
        params (str): Additional parameters required by the executable.
        """

        # Stores the path, executable and parameters.
        self.directory = directory
        self.executable = executable
        self.params = params

    def run(self):
        """
        Run the executable with the given parameters in its directory.
        """

        # Joining the path and the executable name.
        executable_path = os.path.join(self.directory, self.executable)

        # Executing the file on the shell by joining the executable and the parameters.
        subprocess.run([executable_path,*self.params], cwd=self.directory)


