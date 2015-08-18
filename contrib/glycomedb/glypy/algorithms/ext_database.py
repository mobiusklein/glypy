from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import PickleType, Numeric, String, create_engine, Column, Integer, ForeignKey, Text

from glypy import Glycan
from glypy.io import glycoct as glycoct_parser


Base = declarative_base()


class ExtGlycanRecord(Base):
    __tablename__ = "ExtGlycanRecord"
    id = Column(Integer, primary_key=True, autoincrement=True)
    structure = Column(PickleType)
    glycoct = Column(Text, index=True)
    mass = Column(Numeric(12, 4), index=True)
    #taxa = relationship("Taxonomy", backref=backref("glycans", lazy='dynamic'))
    #motifs = relationship("Motif", backref=backref("glycans", lazy='dynamic'))


class Taxon(Base):
    __tablename__ = "Taxonomy"
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column("name")
