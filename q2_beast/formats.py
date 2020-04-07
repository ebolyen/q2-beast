import qiime2.plugin.model as model


class PosteriorLogFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class NexusFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BEASTControlFileFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass

    def md5sum(self):
        import qiime2.core.util  # TODO: don't import from here
        return qiime2.core.util.md5sum(self)


class BEASTOpsFileFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BEASTPosteriorDirFmt(model.DirectoryFormat):
    log = model.File('posterior.log', format=PosteriorLogFormat)
    trees = model.File('posterior.trees', format=NexusFormat)
    ops = model.File('posterior.ops', format=BEASTOpsFileFormat)
    control = model.File('control_file.xml',
                         format=BEASTControlFileFormat)


NexusDirFmt = model.SingleFileDirectoryFormat(
    'NexusDirFmt', 'data.nex', format=NexusFormat)
