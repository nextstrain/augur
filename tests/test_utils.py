import datetime

from augur.utils import ambiguous_date_to_date_range
from augur.utils import run_shell_command

from freezegun import freeze_time
import pytest


@pytest.fixture
def mock_check_output(mocker):
    return mocker.patch("subprocess.check_output")


@pytest.fixture
def mock_check_output_raise_called_process_error(mocker):
    import subprocess

    return mocker.patch(
        "subprocess.check_output", side_effect=subprocess.CalledProcessError(5, 6)
    )


@pytest.fixture
def mock_check_output_raise_file_not_found_error(mocker):
    return mocker.patch("subprocess.check_output", side_effect=FileNotFoundError)


class TestUtils:
    def test_ambiguous_date_to_date_range_not_ambiguous(self):
        assert ambiguous_date_to_date_range("2000-03-29", "%Y-%m-%d") == (
            datetime.date(year=2000, month=3, day=29),
            datetime.date(year=2000, month=3, day=29),
        )

    def test_ambiguous_date_to_date_range_ambiguous_day(self):
        assert ambiguous_date_to_date_range("2000-01-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=1, day=31),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month(self):
        assert ambiguous_date_to_date_range("2000-XX-5", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=5),
            datetime.date(year=2000, month=12, day=5),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month_and_day(self):
        assert ambiguous_date_to_date_range("2000-XX-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=12, day=31),
        )

    @freeze_time("2000-02-20")
    def test_ambiguous_date_to_date_range_current_day_limit(self):
        assert ambiguous_date_to_date_range("2000-02-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=2, day=1),
            datetime.date(year=2000, month=2, day=20),
        )

    def test_run_shell_command(self, mocker, mock_check_output):
        run_shell_command("great-command")

        mock_check_output.assert_called_once_with(
            ["/bin/bash", "-c", "set -euo pipefail; great-command"],
            shell=mocker.ANY,
            stderr=mocker.ANY,
            env=mocker.ANY,
        )

    def test_run_shell_command_not_posix(self, mocker, mock_check_output):
        mocker.patch("os.name", "certainly not posix, that's for sure")

        run_shell_command("great-command")

        mock_check_output.assert_called_once_with(
            ["env", "bash", "-c", "set -euo pipefail; great-command"],
            shell=mocker.ANY,
            stderr=mocker.ANY,
            env=mocker.ANY,
        )

    def test_run_shell_command_extra_env(self, mocker, mock_check_output):
        run_shell_command("great-command", False, {"x": 6})

        assert mock_check_output.call_args[1]["env"]["x"] == 6

    def test_run_shell_command_called_process_error(
        self, mocker, mock_check_output_raise_called_process_error
    ):
        print_error = mocker.patch("augur.utils.print_error")

        run_shell_command("great-command")

        print_error.assert_called_with(
            "None\nshell exited 5 when running: great-command"
        )

    def test_run_shell_command_called_process_error_raise(
        self, mock_check_output_raise_called_process_error
    ):
        import subprocess

        with pytest.raises(subprocess.CalledProcessError):
            run_shell_command("great-command", True)

    def test_run_shell_command_called_process_error_no_raise(
        self, mock_check_output_raise_called_process_error
    ):
        assert run_shell_command("great-command", False) is False

    def test_run_shell_command_file_not_found(
        self, mocker, mock_check_output_raise_file_not_found_error
    ):
        print_error = mocker.patch("augur.utils.print_error")

        run_shell_command("great-command")

        assert "Unable to run shell commands using" in print_error.call_args[0][0]

    def test_run_shell_command_file_not_found_raise(
        self, mock_check_output_raise_file_not_found_error
    ):
        with pytest.raises(FileNotFoundError):
            run_shell_command("great-command", True)

    def test_run_shell_command_file_not_found_no_raise(
        self, mock_check_output_raise_file_not_found_error
    ):
        assert run_shell_command("great-command", False) is False
