import pathlib as pl
from pprint import pformat
from typing import Tuple, Type, Union

from pydantic import Field
from pydantic_settings import BaseSettings, PydanticBaseSettingsSource, TomlConfigSettingsSource

toml_file = pl.Path(__file__).resolve().parent / "config.toml"  # May be updated by the create_config function if the user provides a file path to a config file


class Units(BaseSettings, validate_assignment=True):
    orbital_energy_unit: str = Field("eV", description="The unit of the orbital energies that are read from the adf.rkf file(s)")


class PyFragConfig(BaseSettings, validate_assignment=True):
    units: Units

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: Type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> Tuple[PydanticBaseSettingsSource, ...]:
        return (TomlConfigSettingsSource(settings_cls, toml_file=toml_file), init_settings, env_settings, dotenv_settings, file_secret_settings)

    def __str__(self) -> str:
        formatted = self.model_dump(exclude_defaults=False)
        formatted = pformat(formatted)
        return formatted


def create_config(file_path: Union[pl.Path, str, None]) -> PyFragConfig:
    """Creates the config object from the given file path."""
    global toml_file
    toml_file = pl.Path(file_path) if file_path is not None else pl.Path(__file__).resolve().parent / "config.toml"
    return PyFragConfig()  # type: ignore


pyfrag_config = create_config(None)
