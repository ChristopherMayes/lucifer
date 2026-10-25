"""The package without the executable: the generated models against the schema, the input file
through TOML and back, the converter, and the dotted access.

    python -m pytest python/tests -q
"""
import dataclasses
import tomllib

import pytest

from lucifer import Input, schema
from lucifer import models
from lucifer.flat import to_toml
from lucifer.generate import render


def test_models_match_the_schema_field_by_field():
    """Every table's fields are the schema's keys with the schema's defaults, and the
    committed models.py is what the generator writes from the committed schema."""
    sch = schema.load()
    keys = schema.keys(sch)
    for table, cls in (("", Input), *models.TABLES):
        obj = cls()
        have = {f.name: getattr(obj, f.name) for f in dataclasses.fields(cls)
                if not dataclasses.is_dataclass(getattr(obj, f.name))}
        assert set(have) == set(keys[table]), table
        for k, prop in keys[table].items():
            assert have[k] == prop["default"], f"{table}.{k}"
    assert render(sch) == open(models.__file__).read()


def test_an_input_round_trips_through_toml():
    d = Input(lattice="aramis.bmad")
    d.beam.n_particle = 4096
    d.run.out_root = "s7"
    d.beam.resample.n_slice = 3
    d.field.harmonics = [1, 3]
    d.wake.on = True
    text = d.toml()
    # Only the keys set away from the defaults are written, the lattice always.
    assert text.splitlines()[0] == 'lattice = "aramis.bmad"'
    assert "[beam.resample]" in text and "n_slice = 3" in text
    assert "source_filter" not in text
    assert Input.loads(text) == d
    # The explicit form carries every key and reads back to the same input file.
    full = d.toml(explicit=True)
    assert Input.loads(full) == d
    assert len(tomllib.loads(full)["run"]) == len(dataclasses.fields(models.Run))


def test_the_explicit_input_names_every_key_once():
    data = tomllib.loads(Input(lattice="a.bmad").toml(explicit=True))
    keys = schema.keys()
    for table, names in keys.items():
        tab = data
        for part in [p for p in table.split(".") if p]:
            tab = tab[part]
        own = {k for k, v in tab.items() if not isinstance(v, dict)}
        assert own == set(names), table


def test_from_flat_agrees_with_the_converter():
    flat = "lat_file = 'a.bmad'\nglobal%out_root = 'q'\nn_particle = 2048\nlambda0 = 1e-10\nharmonics = 1, 3\n"
    d = Input.from_flat(flat)
    assert d == Input.loads(to_toml(flat))
    assert d.lattice == "a.bmad" and d.run.out_root == "q"
    assert d.beam.n_particle == 2048 and d.field.lambda0 == 1e-10 and d.field.harmonics == [1, 3]


def test_dotted_access_reaches_the_field():
    d = Input()
    d["beam.n_particle"] = 512
    d["beam.resample.slice_width"] = 0.02
    d["lattice"] = "x.bmad"
    assert d.beam.n_particle == 512 and d.beam.resample.slice_width == 0.02 and d.lattice == "x.bmad"
    assert d["beam.resample.slice_width"] == 0.02
    assert ("beam.n_particle", 512) in d.items()
    with pytest.raises(KeyError):
        d["beam.no_such_key"]


def test_an_unknown_key_is_pythons_own_error():
    with pytest.raises(TypeError):
        Input.loads('lattice = "a"\n[beam]\nlambda0 = 1e-10\n')
    with pytest.raises(TypeError):
        models.Run(no_such=1)
