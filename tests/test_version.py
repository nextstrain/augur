from augur import version


class MockVersion:

    def get_version(self):
        return self.version

    def set_version(self, version):
        self.version = version

    def print_version(self):
        print(self.version)


class TestVersion:

    def test_get_version(self):
        MockVersion.set_version(self, '0.1.0')
        self.mock_version = MockVersion.get_version(self)
        assert self.mock_version.split('.')[0] == '0'
        assert self.mock_version.split('.')[1] == '1'
        assert self.mock_version.split('.')[2] == '0'

    def test_get_version_is_not_none(self):
        assert version._get_version() is not None

    def test_run_with_no_args(self):
        args = None
        assert version.run(args) == 0
