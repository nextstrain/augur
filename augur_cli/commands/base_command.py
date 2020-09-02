from abc import ABC, abstractmethod
import re


class BaseCommand(ABC):
    def __init__(self, args):
        self.args = args

    @abstractmethod
    def run(self):
        pass

    #TODO other abstractmethods

    @classmethod
    def command_name(cls):
        cls_name_without_command = cls.__name__[0:-7]
        return re.sub("(?!^)([A-Z])", r"-\1", cls_name_without_command).lower()
