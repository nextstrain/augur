import pytest
import augur.io.file


class TestFile:
    def test_open_file_read_text(self, tmpdir):
        """Read a text file."""
        path = str(tmpdir / 'test.txt')
        with open(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with augur.io.file.open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_text(self, tmpdir):
        """Write a text file."""
        path = str(tmpdir / 'test.txt')
        with augur.io.file.open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with open(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_gzip(self, tmpdir):
        """Read a text file compressed with gzip."""
        import gzip
        path = str(tmpdir / 'test.txt.gz')
        with gzip.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')
        with augur.io.file.open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_gzip(self, tmpdir):
        """Write a text file compressed with gzip."""
        import gzip
        path = str(tmpdir / 'test.txt.gz')
        with augur.io.file.open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with gzip.open(path, 'rt') as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_lzma(self, tmpdir):
        """Read a text file compressed with LZMA."""
        import lzma
        path = str(tmpdir / 'test.txt.xz')
        with lzma.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')
        with augur.io.file.open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_lzma(self, tmpdir):
        """Write a text file compressed with LZMA."""
        import lzma
        path = str(tmpdir / 'test.txt.xz')
        with augur.io.file.open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with lzma.open(path, 'rt') as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_zstd(self, tmpdir):
        """Read a text file compressed with zstd."""
        import zstandard as zstd
        path = str(tmpdir / 'test.txt.zst')
        with zstd.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')
        with augur.io.file.open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_zstd(self, tmpdir):
        """Write a text file compressed with zstd."""
        import zstandard as zstd
        path = str(tmpdir / 'test.txt.zst')
        with augur.io.file.open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with zstd.open(path, 'rt') as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_nested_read(self, tmpdir):
        """Open a text file as a path then buffer for reading."""
        path = str(tmpdir / 'test.txt')
        with open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')

        with augur.io.file.open_file(path) as f_read:
            with augur.io.file.open_file(f_read) as f_read:
                assert f_read.read() == 'foo\nbar\n'

    def test_open_file_nested_read_zstd(self, tmpdir):
        """Open a zstd-compressed text file as a path then buffer for reading."""
        import zstandard as zstd
        path = str(tmpdir / 'test.txt.zst')
        with zstd.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')

        with augur.io.file.open_file(path) as f_read:
            with augur.io.file.open_file(f_read) as f_read:
                assert f_read.read() == 'foo\nbar\n'

    def test_open_file_nested_write(self, tmpdir):
        """Open a text file as a path then buffer for writing."""
        path = str(tmpdir / 'test.txt')
        with augur.io.file.open_file(path, 'w') as f_write:
            with augur.io.file.open_file(f_write, 'w') as f_write:
                f_write.write('foo\nbar\n')

    def test_open_file_nested_write_zstd(self, tmpdir):
        """Open a zstd-compressed text file as a path then buffer for writing."""
        path = str(tmpdir / 'test.txt.zst')
        with augur.io.file.open_file(path, 'w') as f_write:
            with augur.io.file.open_file(f_write, 'w') as f_write:
                f_write.write('foo\nbar\n')

    def test_open_file_read_error(self):
        """Try reading from an unsupported type."""
        path_or_buffer = len("bogus")

        with pytest.raises(TypeError, match="Type <class 'int'> is not supported."):
            with augur.io.file.open_file(path_or_buffer, 'r') as f:
                f.read()
