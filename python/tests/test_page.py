"""The examples on docs/python.md that need no executable, run as tests.

Each region between "# example: <name>" and "# end example: <name>" is the code of one
block on the page, byte for byte after the common indentation is removed, which
tests/scripts/check_python_page.py holds. The assertions after a region say what the
page's prose says the example does. The examples that run the program are in the
package section of tests/scripts/check_program.py, which the keystone runs.
"""
import tomllib

from lucifer import Input


def test_input(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # example: input
    from lucifer import Input

    inp = Input(lattice="aramis.bmad")
    inp.beam.n_particle = 4096
    inp["run.out_root"] = "s7"          # a key by its dotted name
    inp.write("s7.toml")                # the TOML the program reads
    same = Input.read("s7.toml")
    # end example: input
    assert same == inp
    # Only the keys set away from the defaults are written, the lattice always.
    assert tomllib.loads((tmp_path / "s7.toml").read_text()) == {
        "lattice": "aramis.bmad", "run": {"out_root": "s7"}, "beam": {"n_particle": 4096}}


def test_explicit(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    inp = Input(lattice="aramis.bmad")
    # example: explicit
    inp.write("s7-all.toml", explicit=True)   # every key, the defaults written out
    for key, value in inp.items():            # the keys set away from their defaults
        print(key, "=", value)
    # end example: explicit
    data = tomllib.loads((tmp_path / "s7-all.toml").read_text())
    assert data["run"]["ran_seed"] == 12345 and inp.items() == [("lattice", "aramis.bmad")]


def test_flat():
    # example: flat
    from lucifer.flat import to_toml

    text = to_toml("lattice = 'aramis.bmad'\nn_particle = 2048\nrun.out_root = 's7'\n")
    # end example: flat
    assert tomllib.loads(text) == {"lattice": "aramis.bmad", "run": {"out_root": "s7"},
                                   "beam": {"n_particle": 2048}}


def test_schema():
    # example: schema
    from lucifer import schema

    schema.path()                                   # schema.json, as the program wrote it
    schema.keys()["run"]["ran_seed"]["default"]     # 12345
    # end example: schema
    assert schema.path().name == "schema.json"
    assert schema.keys()["run"]["ran_seed"]["default"] == 12345
