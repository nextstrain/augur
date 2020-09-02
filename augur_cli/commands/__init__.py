import importlib
import pkgutil

from augur_cli.commands.base_command import BaseCommand

for module_info in pkgutil.iter_modules(__path__, prefix=f"{__name__}."):
    importlib.import_module(module_info.name)

COMMAND_CLASSES = [command_class for command_class in BaseCommand.__subclasses__()]
