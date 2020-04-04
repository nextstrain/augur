import factory
import inspect
import pytest_factoryboy

import tests.factories

# register all of the factory.Factory subclasses in tests.factories
for (name, factory_class) in inspect.getmembers(tests.factories):
    if inspect.isclass(factory_class) and issubclass(factory_class, factory.Factory):
        pytest_factoryboy.register(factory_class)
