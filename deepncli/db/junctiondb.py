from peewee import *
from playhouse.apsw_ext import APSWDatabase
from playhouse.migrate import SqliteMigrator, migrate


database = APSWDatabase(None, pragmas={'journal_mode': 'off',
                                       'cache_size': -500*1000,
                                       'ignore_check_constraints': 0,
                                       'synchronous': 0,
                                       'foreign_keys': 1,
                                       'temp_store': 'memory',
                                       'locking_mode': 'exclusive'})  # Un-initialized database.


class Gene(Model):
    gene_name = TextField()
    orf_start = IntegerField()
    orf_stop = IntegerField()
    mrna = TextField()
    intron = CharField()
    chromosome = CharField()
    nm_number = CharField()

    class Meta:
        database = database


class Junction(Model):
    gene = ForeignKeyField(Gene)
    position = IntegerField()
    query_start = IntegerField()
    frame = TextField()
    ppm = FloatField()
    orf = TextField()
    inframe_inorf = BooleanField()
    count = IntegerField()

    class Meta:
        database = database


class Stats(Model):
    gene = ForeignKeyField(Gene)
    backwards = IntegerField(default=0)
    downstream = IntegerField(default=0)
    inframe_inorf = IntegerField(default=0)
    in_frame = IntegerField(default=0)
    in_orf = IntegerField(default=0)
    intron = IntegerField(default=0)
    not_in_frame = IntegerField(default=0)
    total = IntegerField(default=0)
    upstream = IntegerField(default=0)

    class Meta:
        database = database


class JunctionsDatabase(object):
    def __init__(self, db_name):
        self.db = database
        self.db.init(db_name)
        self.migrator = SqliteMigrator(self.db)

    def create_tables(self):
        self.db.drop_tables([Stats, Junction, Gene])
        self.db.create_tables([Gene, Junction, Stats])
        self.create_indexes()

    def create_indexes(self):
        migrate(self.migrator.add_index('gene', ('id', 'gene_name', 'nm_number')),
                self.migrator.add_index('junction', ('id', 'gene_id')),
                self.migrator.add_index('stats', ('id', 'gene_id')))

    def close_db(self):
        self.db.close()
