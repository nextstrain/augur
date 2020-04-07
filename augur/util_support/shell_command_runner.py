import os
import sys
import subprocess
from textwrap import dedent


class ShellCommandRunner:
    """
    Run the given command string via Bash with error checking.

    TODO move to method docstrings
    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.

    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """

    def __init__(self, cmd, *, raise_errors=False, extra_env=None):
        self.cmd = cmd
        self.raise_errors = raise_errors
        self.extra_env = extra_env

    def run(self):
        try:
            self.invoke_command()
        except Exception as error:
            self.print_error_message(error)

            if self.raise_errors:
                raise error

            return False

        return True

    def invoke_command(self):
        return subprocess.check_output(
            self.shell_executable + self.shell_args,
            shell=False,
            stderr=subprocess.STDOUT,
            env=self.modified_env,
        )

    @property
    def shell_executable(self):
        if os.name == "posix":
            return ["/bin/bash"]
        else:
            # We try best effort on other systems. For now that means nt/java.
            return ["env", "bash"]

    @property
    def shell_args(self):
        return ["-c", "set -euo pipefail; " + self.cmd]

    @property
    def modified_env(self):
        env = os.environ.copy()

        if self.extra_env:
            env.update(self.extra_env)

        return env

    def print_error_message(self, error):
        if isinstance(error, subprocess.CalledProcessError):
            message = f"{error.output}\nshell exited {error.returncode} when running: {self.cmd}"

            if error.returncode == 127:
                message += "\nAre you sure this program is installed?"
        elif isinstance(error, FileNotFoundError):
            shell = " and ".join(self.shell_executable)

            message = f"""
                Unable to run shell commands using {shell}!

                Augur requires {shell} to be installed.  Please open an issue on GitHub
                <https://github.com/nextstrain/augur/issues/new> if you need assistance.
                """
        else:
            message = str(error)

        self.print_error(message)

    @staticmethod
    def print_error(message):
        """Prints message to STDERR formatted with textwrap.dedent"""
        print("\nERROR: " + dedent(message).lstrip("\n") + "\n", file=sys.stderr)
