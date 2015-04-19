from sqlalchemy.ext.declarative import declarative_base, DeclarativeMeta
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import PickleType, Numeric, String, create_engine, Column, Integer, ForeignKey, Text

from pygly2 import Glycan
from pygly2.io import glycoct as glycoct_parser

Base = declarative_base()

class ExtGlycanRecord(Base):
    __tablename__ = "ExtGlycanRecord"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    structure = Column("structure", PickleType)
    glycoct = Column("glycoct", Text)
    mass = Column("mass", Numeric)
    taxa = relationship("Taxonomy", backref=backref("glycans"))
    motifs = relationship("Motif", backref=backref("glycans"))


class Taxon(Base):
    __tablename__ = "Taxonomy"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
